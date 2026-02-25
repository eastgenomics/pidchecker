#!/usr/bin/env python3
"""
PII Checker for Genomic Variant Classification Comments
========================================================
Uses Microsoft Presidio with custom recognizers tuned for the UK clinical
genomics context.  Produces an HTML report and a CSV summary of findings.

Setup
-----
    pip install -r requirements.txt
    python -m spacy download en_core_web_lg
    python pii_checker.py
"""

import csv
import html as html_lib
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

from presidio_analyzer import (
    AnalyzerEngine,
    Pattern,
    PatternRecognizer,
    RecognizerRegistry,
    RecognizerResult,
)
from presidio_analyzer.nlp_engine import NlpEngineProvider


# ──────────────────────────────────────────────────────────────────────────────
# Configuration  (edit these to suit your environment)
# ──────────────────────────────────────────────────────────────────────────────

INPUT_FILE   = Path("comment.csv")
OUTPUT_HTML  = Path("pii_report.html")
OUTPUT_CSV   = Path("pii_findings.csv")

# Minimum Presidio confidence score to include in output (0.0 - 1.0).
# Lower values catch more but produce more false positives.
SCORE_THRESHOLD = 0.35

# spaCy model.  en_core_web_lg gives better NER than _sm; _trf is best but slow.
SPACY_MODEL = "en_core_web_lg"

# Entities to scan for.  Custom entities (UK_NHS etc.) are added below.
# DATE_TIME and LOCATION are intentionally excluded:
#   DATE_TIME - generates too many false positives from publication dates and
#               database timestamps; a dedicated DOB recognizer is used instead.
#   LOCATION  - "European", "South Asian" etc. are population descriptors, not PII.
ENTITIES = [
    "PERSON",
    "EMAIL_ADDRESS",
    "URL",
    "UK_NHS",         # custom - NHS number with modulus-11 checksum
    "UK_POSTCODE",    # custom
    "DATE_OF_BIRTH",  # custom - dates preceded by DOB context words only
    "PATIENT_ID",     # custom - adjust regex to match your internal ID format
]

# Colours used to highlight entity types in the HTML report
ENTITY_COLOURS: Dict[str, str] = {
    "PERSON":        "#ff9999",
    "EMAIL_ADDRESS": "#ffcc80",
    "URL":           "#e0e0e0",
    "UK_NHS":        "#ff6680",
    "UK_POSTCODE":   "#b39ddb",
    "DATE_OF_BIRTH": "#f48fb1",
    "PATIENT_ID":    "#ff8f00",
}
_DEFAULT_COLOUR = "#ffe082"


# ──────────────────────────────────────────────────────────────────────────────
# NHS Number Validation
# ──────────────────────────────────────────────────────────────────────────────

def _validate_nhs(text: str) -> bool:
    """Return True if *text* passes the NHS number modulus-11 check digit."""
    digits = re.sub(r"\D", "", text)
    if len(digits) != 10:
        return False
    total = sum(int(d) * (10 - i) for i, d in enumerate(digits[:9]))
    remainder = 11 - (total % 11)
    if remainder == 11:
        check_digit = 0
    elif remainder == 10:
        return False          # No valid NHS number has this remainder
    else:
        check_digit = remainder
    return int(digits[9]) == check_digit


# ──────────────────────────────────────────────────────────────────────────────
# Custom Recognizers
# ──────────────────────────────────────────────────────────────────────────────

def _nhs_recognizer() -> PatternRecognizer:
    """
    Match 10-digit NHS numbers in both spaced (XXX XXX XXXX) and compact forms.
    Checksum validation is applied during post-processing filtering.
    The compact form receives a lower initial score because bare 10-digit
    sequences are more ambiguous (e.g. allele counts, ClinVar IDs).
    """
    return PatternRecognizer(
        supported_entity="UK_NHS",
        patterns=[
            Pattern("NHS_SPACED",  r"\b\d{3} \d{3} \d{4}\b", 0.85),
            Pattern("NHS_COMPACT", r"\b\d{10}\b",             0.50),
        ],
    )


def _postcode_recognizer() -> PatternRecognizer:
    return PatternRecognizer(
        supported_entity="UK_POSTCODE",
        patterns=[
            Pattern(
                "UK_POSTCODE",
                r"(?i)\b[A-Z]{1,2}[0-9][0-9A-Z]?\s?[0-9][A-Z]{2}\b",
                0.75,
            ),
        ],
    )


def _dob_recognizer() -> PatternRecognizer:
    """
    Flag dates only when they appear alongside explicit DOB context words.
    This avoids the large number of false positives from publication years
    and database timestamps (e.g. "as of 2023-10-05").
    The bare DD/MM/YYYY pattern is retained at a low score so that
    uncontextualised dates can still surface for manual review.
    """
    return PatternRecognizer(
        supported_entity="DATE_OF_BIRTH",
        patterns=[
            Pattern(
                "DOB_LABELLED",
                r"(?i)(?:DOB|D\.O\.B\.?|date\s+of\s+birth|born|dob\s*:)\s*"
                r"\d{1,2}[\/\-\.]\d{1,2}[\/\-\.]\d{2,4}",
                0.95,
            ),
            Pattern(
                "DATE_DD_MM_YYYY",
                r"\b\d{1,2}[\/\-]\d{1,2}[\/\-]\d{4}\b",
                0.40,
            ),
        ],
    )


def _patient_id_recognizer() -> PatternRecognizer:
    """
    Recognizer for internal patient / sample IDs.
    Add further Pattern entries here as new ID formats are identified.
    """
    return PatternRecognizer(
        supported_entity="PATIENT_ID",
        patterns=[
            Pattern(
                "PATIENT_ID_PREFIXED",
                r"\b(?:PAT|LAB|SMP|REF|SAMPLE|PATIENT)[-\s]?\d{4,10}\b",
                0.85,
            ),
            Pattern(
                "SPECIMEN_ACCESSION",   # e.g. SP12345R0123
                r"\bSP\d+R\d+\b",
                0.90,
            ),
            Pattern(
                "GM_CASE_NUMBER",       # e.g. GM1234567
                r"\bGM\d+\b",
                0.90,
            ),
        ],
    )


# ──────────────────────────────────────────────────────────────────────────────
# False-Positive Filtering
# ──────────────────────────────────────────────────────────────────────────────

# ── Safe-span patterns ────────────────────────────────────────────────────────
# Any Presidio finding whose span overlaps with a match from one of these
# patterns is suppressed as a known non-PII false positive.

# ACMG/AMP variant classification criteria codes (standalone or with strength).
# Matches: PM2, PVS1_Strong, PS3 - Strong, PP4 - Supporting, BA1 - Standalone …
_ACMG_CODE = re.compile(
    r"\b(?:PVS|PS|PM|PP|BS|BP|BA)\d+(?:_\w+)?"
    r"(?:\s*[-\u2013]\s*(?:Very\s+Strong|Strong|Moderate|Supporting|Standalone))?\b",
    re.IGNORECASE,
)

# HGVS variant nomenclature and RefSeq accession numbers.
# Matches: c.1234A>G, p.Arg3527Gln, NM_000527.4, NM_000527.4(LDLR):c.1730G>C,
#          NP_005903.2:p.Y35*, p700Asn (non-standard), c.*60_*61insCTTTA, etc.
_HGVS_NOTATION = re.compile(
    r"(?:"
    # RefSeq accessions with optional gene symbol and HGVS suffix
    r"\b(?:NM|NP|NR|NG|NC|NT|LRG|XM|XP)_\d+(?:\.\d+)?(?:\([A-Z0-9]+\))?"
    r"(?::[cgnrmp]\.[^\s,;\[\]()]+)?"
    # cDNA / RNA / genomic: c.1234A>G, c.*60_*61ins…, g.12345A>G
    r"|(?<![A-Za-z])[cgnrm]\.[*+\-0-9][^\s,;\[\]]*"
    # Standard protein: p.Arg3527Gln, p.Val63Ile, p.Y35*, p.Met1?
    r"|(?<![A-Za-z])p\.[A-Z*][^\s,;\[\]]*"
    # Non-standard protein without dot: p700Asn
    r"|\bp\d+[A-Z][a-z]{2}\b"
    r")",
    re.UNICODE,
)

# Literature citations: "Surname et al. YEAR, Journal Vol:Page"
# Uses the Latin Extended Unicode block (U+00C0–U+024F) so accented surnames
# such as Bañares, Graça, Böhm, etc. are all matched.
_CITATION_REF = re.compile(
    r"\b[A-Z\u00C0-\u024F][a-z\u00C0-\u024F\-]+"        # first surname
    r"(?:\s+[A-Z\u00C0-\u024F][a-z\u00C0-\u024F\-]+)*"  # optional further name parts
    r"\s+et\s+al\."                                       # "et al."
    r"(?:\s+\d{4}(?:,\s+[^;)\n]+)?)?",                   # optional: year + journal ref
    re.UNICODE,
)


# Genomic variant database and population resource names.
# gnomAD allele-count strings like "1/113.664 (0.0009%)" are also handled
# because gnomAD itself is in this list and the surrounding numeric context
# will not match any PII entity.
_DATABASE_NAME = re.compile(
    r"\b(?:gnomAD|ClinVar|1kG|dbSNP|ExAC|HGMD|LOVD|OMIM|UniProt|dbVar|"
    r"1000\s+Genomes|Ensembl|RefSeq|GnomAD)\b",
    re.IGNORECASE,
)

# HGNC-approved gene symbols commonly encountered in clinical genetics reporting.
# Extend this list as required by your panel scope.
_GENE_SYMBOL = re.compile(
    r"\b(?:"
    r"ABCA4|ABCB11|ABCC8|ABCD1|AGL|AHI1|ALDOB|ALPL|APC|APOB|APOC2|APOE|AR|ARX|"
    r"ASS1|ATM|ATP7A|ATP7B|"
    r"BBS1|BCOR|BEST1|BRCA1|BRCA2|"
    r"CACNA1A|CASR|CDH1|CFTR|CHD7|CHRNE|COL1A1|COL1A2|COL4A3|COL4A4|COL4A5|CPS1|"
    r"CYP27A1|CYP2R1|CYP7B1|"
    r"DHCR7|DMD|DNMT3B|DUOX2|"
    r"EDA|EFNB1|ELN|EMD|EPM2A|EXT1|EXT2|"
    r"FAH|FANCA|FANCC|FGFR3|FMR1|FOXG1|FOXP3|"
    r"G6PC|GALT|GBA|GDF5|GJB1|GJB2|GJB6|GLA|GNAS|"
    r"HBA1|HBA2|HBB|HFE|HNF1A|HNF1B|HNF4A|HRAS|"
    r"IDS|IDUA|IKBKG|IL2RG|INSR|"
    r"JAG1|"
    r"KAL1|KCNJ11|KCNQ1|KCNQ4|KDM5C|KRT1|KRT10|KRT14|"
    r"L1CAM|LAMP2|LDLR|LEPR|LIPA|LPL|"
    r"MC4R|MECP2|MEN1|MLH1|MPZ|MSH2|MSH6|MYT1L|"
    r"NDP|NF1|NF2|NOTCH3|"
    r"OCA2|OTC|"
    r"PAH|PCDH19|PCSK9|PDHA1|PHEX|PHIP|PKD1|PKD2|PKP2|PMP22|POMT1|POMT2|"
    r"PRKAR1A|PRKCG|PRSS1|PTCH1|PTEN|"
    r"RAB7A|RB1|RET|RPE65|RYR1|"
    r"SCN1A|SCN2A|SCN5A|SCN9A|SETD5|SHANK3|SLC25A13|SLC2A1|SLC7A9|SMAD3|SMAD4|"
    r"SOD1|SPINK1|STAT3|STK11|"
    r"TAZ|TBX5|TGFBR1|TGFBR2|TMPRSS3|TP53|TREX1|TSC1|TSC2|"
    r"UBE3A|UGT1A1|"
    r"VHL|VPS13A|"
    r"WFS1|WT1|"
    r"ZEB2"
    r")\b",
    re.UNICODE,
)

# Common biomedical journal abbreviations used in clinical genetics literature.
# Matched with an optional trailing volume:page (e.g. "J Med Genet 40:123").
_JOURNAL_NAME = re.compile(
    r"\b(?:"
    r"Am\s+J\s+(?:Hum\s+)?Genet|Am\s+J\s+Med\s+Genet(?:\s+[A-C])?|"
    r"Ann\s+(?:Hum\s+)?Genet|Ann\s+Intern\s+Med|Ann\s+Neurol|"
    r"Atherosclerosis|Biomedicines|Brain|"
    r"Circ\s+Genom\s+Precis\s+Med|Circ\s+Genet|"
    r"Clin\s+(?:Biochem|Chem|Genet)|"
    r"Eur\s+Heart\s+J|Eur\s+J\s+Hum\s+Genet|"
    r"Genet\s+(?:Med|Test)|"
    r"Heart|Hum\s+(?:Genet|Mol\s+Genet|Mutat|Reprod)|"
    r"J\s+Am\s+Coll\s+Cardiol|"
    r"J\s+Clin\s+(?:Endocrinol\s+Metab|Lipidol|Invest)|"
    r"J\s+Genet\s+Couns|J\s+Hum\s+Genet|J\s+Inherit\s+Metab\s+Dis|"
    r"J\s+Med\s+Genet|J\s+Neurol(?:\s+Sci)?|J\s+Pediatr|"
    r"Lancet|"
    r"Mol\s+Genet\s+(?:Genomic\s+Med|Metab)|"
    r"N\s+Engl\s+J\s+Med|Nat\s+(?:Genet|Rev\s+Genet|Med)|"
    r"Neurology|Neuromuscul\s+Disord|NPJ\s+Genom\s+Med|"
    r"Orphanet\s+J\s+Rare\s+Dis|"
    r"Pediatrics|PLOS\s+(?:Genet|One)|"
    r"Sci\s+Rep|Stroke"
    r")(?:\s+\d+[:\d]*)?",
    re.IGNORECASE,
)


def _safe_spans(text: str) -> List[Tuple[int, int]]:
    """Return (start, end) spans of text regions known to be non-PII."""
    spans: List[Tuple[int, int]] = []
    for pattern in (
        _ACMG_CODE,
        _HGVS_NOTATION,
        _CITATION_REF,
        _DATABASE_NAME,
        _GENE_SYMBOL,
        _JOURNAL_NAME,
    ):
        for m in pattern.finditer(text):
            spans.append((m.start(), m.end()))
    return spans


def _overlaps_any(start: int, end: int, safe_spans: List[Tuple[int, int]]) -> bool:
    """Return True if [start, end) overlaps with any safe span."""
    return any(start < se and end > ss for ss, se in safe_spans)


def _filter_results(
    results: List[RecognizerResult], text: str
) -> List[RecognizerResult]:
    """Remove known false positives from the raw Presidio output."""
    safe = _safe_spans(text)
    kept: List[RecognizerResult] = []
    for r in results:
        # UK_NHS: reject if the match fails the modulus-11 checksum
        if r.entity_type == "UK_NHS":
            if not _validate_nhs(text[r.start: r.end]):
                continue

        # Suppress anything overlapping a known-safe span
        # (ACMG criteria codes, HGVS nomenclature, literature citations)
        if _overlaps_any(r.start, r.end, safe):
            continue

        kept.append(r)
    return kept


# ──────────────────────────────────────────────────────────────────────────────
# Analyzer Construction
# ──────────────────────────────────────────────────────────────────────────────

def build_analyzer() -> AnalyzerEngine:
    """Assemble the Presidio AnalyzerEngine with all recognizers."""
    nlp_config = {
        "nlp_engine_name": "spacy",
        "models": [{"lang_code": "en", "model_name": SPACY_MODEL}],
    }
    try:
        nlp_engine = NlpEngineProvider(nlp_configuration=nlp_config).create_engine()
    except OSError:
        print(
            f"\nERROR: spaCy model '{SPACY_MODEL}' not found.\n"
            f"Install it with:  python -m spacy download {SPACY_MODEL}\n",
            file=sys.stderr,
        )
        sys.exit(1)

    registry = RecognizerRegistry()
    registry.load_predefined_recognizers()
    for rec in [
        _nhs_recognizer(),
        _postcode_recognizer(),
        _dob_recognizer(),
        _patient_id_recognizer(),
    ]:
        registry.add_recognizer(rec)

    return AnalyzerEngine(nlp_engine=nlp_engine, registry=registry)


# ──────────────────────────────────────────────────────────────────────────────
# Per-comment Analysis
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class Finding:
    row: int
    entity_type: str
    matched_text: str
    score: float
    start: int
    end: int
    snippet: str          # surrounding context for the HTML / CSV report


def analyze_comment(
    analyzer: AnalyzerEngine, text: str, row: int
) -> List[Finding]:
    if not text.strip():
        return []

    raw: List[RecognizerResult] = analyzer.analyze(
        text=text,
        entities=ENTITIES,
        language="en",
        score_threshold=SCORE_THRESHOLD,
    )
    raw = _filter_results(raw, text)

    findings: List[Finding] = []
    for r in raw:
        s = max(0, r.start - 60)
        e = min(len(text), r.end + 60)
        snippet = (
            ("…" if s > 0 else "")
            + text[s:e]
            + ("…" if e < len(text) else "")
        )
        findings.append(
            Finding(
                row=row,
                entity_type=r.entity_type,
                matched_text=text[r.start: r.end],
                score=r.score,
                start=r.start,
                end=r.end,
                snippet=snippet,
            )
        )
    return findings


# ──────────────────────────────────────────────────────────────────────────────
# HTML Report
# ──────────────────────────────────────────────────────────────────────────────

def _highlight_text(text: str, findings: List[Finding]) -> str:
    """
    Return the comment as HTML with PII spans wrapped in <mark> elements.
    Overlapping spans are resolved by taking the highest-scoring match.
    """
    if not findings:
        return html_lib.escape(text)

    # Resolve overlaps by score (highest confidence wins), then sort for rendering.
    by_score = sorted(findings, key=lambda f: (-f.score, f.start, f.end))
    chosen: List[Finding] = []
    for f in by_score:
        if not any(f.start < c.end and f.end > c.start for c in chosen):
            chosen.append(f)
    non_overlapping = sorted(chosen, key=lambda f: f.start)

    parts: List[str] = []
    cursor = 0
    for f in non_overlapping:
        parts.append(html_lib.escape(text[cursor: f.start]))
        colour = ENTITY_COLOURS.get(f.entity_type, _DEFAULT_COLOUR)
        label  = html_lib.escape(f.entity_type)
        match  = html_lib.escape(text[f.start: f.end])
        parts.append(
            f'<mark style="background:{colour};border-radius:3px;padding:1px 4px;" '
            f'title="{label} — score: {f.score:.2f}">'
            f"{match}"
            f'<sup style="font-size:0.65em;color:#444;margin-left:2px">'
            f"{label}</sup>"
            f"</mark>"
        )
        cursor = f.end
    parts.append(html_lib.escape(text[cursor:]))
    return "".join(parts)


def _badge(entity_type: str, score: float) -> str:
    colour = ENTITY_COLOURS.get(entity_type, _DEFAULT_COLOUR)
    return (
        f'<span style="background:{colour};padding:1px 7px;border-radius:3px;'
        f'font-size:0.8em;white-space:nowrap">'
        f"{html_lib.escape(entity_type)} {score:.2f}</span>"
    )


def generate_html(
    rows: List[Tuple[int, str, List[Finding]]],
    all_findings: List[Finding],
    total_comments: int,
    output_path: Path,
) -> None:
    flagged_count = sum(1 for _, _, f in rows if f)

    # Entity-type count table
    entity_counts: Dict[str, int] = {}
    for f in all_findings:
        entity_counts[f.entity_type] = entity_counts.get(f.entity_type, 0) + 1
    summary_rows = "".join(
        f"<tr><td>{html_lib.escape(e)}</td><td>{n}</td></tr>"
        for e, n in sorted(entity_counts.items())
    )

    # Legend
    legend = " ".join(
        f'<span style="background:{c};padding:2px 10px;border-radius:3px;'
        f'margin:2px;display:inline-block">{html_lib.escape(e)}</span>'
        for e, c in ENTITY_COLOURS.items()
    )

    # One card per flagged comment
    cards: List[str] = []
    for row_num, text, findings in rows:
        if not findings:
            continue
        badges = " ".join(_badge(f.entity_type, f.score) for f in findings)
        highlighted = _highlight_text(text, findings)
        cards.append(
            f"""
            <div class="card">
              <div class="card-header">Row {row_num}&ensp;{badges}</div>
              <div class="card-body">{highlighted}</div>
            </div>"""
        )

    cards_html = "\n".join(cards) if cards else "<p>No findings.</p>"

    html_out = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>PII Check Report</title>
  <style>
    body  {{ font-family: system-ui, sans-serif; max-width: 1100px;
             margin: 2em auto; padding: 0 1em; color: #222; }}
    h1    {{ color: #b71c1c; }}
    h2    {{ color: #444; margin-top: 2em; border-bottom: 1px solid #ddd;
             padding-bottom: 0.3em; }}
    .grid {{ display: grid; grid-template-columns: repeat(3,1fr);
             gap: 1em; margin: 1em 0; }}
    .stat {{ background: #f5f5f5; border-radius: 6px; padding: 1em;
             text-align: center; }}
    .num  {{ font-size: 2.2em; font-weight: bold; color: #b71c1c; }}
    .lbl  {{ font-size: 0.85em; color: #666; }}
    table {{ border-collapse: collapse; width: auto; }}
    th, td {{ border: 1px solid #ddd; padding: 6px 14px; text-align: left; }}
    th    {{ background: #f5f5f5; }}
    .card {{ border: 1px solid #ddd; border-radius: 6px; margin: 1em 0;
             overflow: hidden; }}
    .card-header {{ background: #fafafa; padding: 0.5em 1em;
                    font-size: 0.85em; border-bottom: 1px solid #ddd; }}
    .card-body {{ padding: 0.8em 1em; font-size: 0.88em; line-height: 1.75;
                  word-break: break-word; }}
  </style>
</head>
<body>
  <h1>PII Check Report — Genomic Variant Comments</h1>

  <h2>Summary</h2>
  <div class="grid">
    <div class="stat"><div class="num">{total_comments}</div>
      <div class="lbl">Comments checked</div></div>
    <div class="stat"><div class="num">{flagged_count}</div>
      <div class="lbl">Comments flagged</div></div>
    <div class="stat"><div class="num">{len(all_findings)}</div>
      <div class="lbl">Total findings</div></div>
  </div>

  <h2>Findings by Entity Type</h2>
  <table>
    <tr><th>Entity type</th><th>Count</th></tr>
    {summary_rows if summary_rows else "<tr><td colspan='2'>None</td></tr>"}
  </table>

  <h2>Legend</h2>
  <p>{legend}</p>
  <p style="font-size:0.85em;color:#666">
    Hover over a highlighted span to see the entity type and confidence score.
    A score of 1.0 is highest confidence. Findings below
    {SCORE_THRESHOLD} are suppressed.
  </p>

  <h2>Flagged Comments</h2>
  {cards_html}
</body>
</html>
"""
    output_path.write_text(html_out, encoding="utf-8")
    print(f"  HTML report  → {output_path}")


# ──────────────────────────────────────────────────────────────────────────────
# CSV Summary
# ──────────────────────────────────────────────────────────────────────────────

def write_csv(findings: List[Finding], output_path: Path) -> None:
    with output_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["row", "entity_type", "score", "matched_text", "snippet"],
        )
        writer.writeheader()
        for f in findings:
            writer.writerow(
                {
                    "row": f.row,
                    "entity_type": f.entity_type,
                    "score": f"{f.score:.3f}",
                    "matched_text": f.matched_text,
                    "snippet": f.snippet,
                }
            )
    print(f"  CSV findings → {output_path}")


# ──────────────────────────────────────────────────────────────────────────────
# Entry Point
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    print(f"Loading spaCy model '{SPACY_MODEL}' and building analyzer…")
    analyzer = build_analyzer()
    print("Analyzer ready.\n")

    # ── Read comments ─────────────────────────────────────────────────────────
    comments: List[Tuple[int, str]] = []
    required_column = "comment_on_classification"
    with INPUT_FILE.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames or required_column not in reader.fieldnames:
            print(
                f"ERROR: required column '{required_column}' not found in {INPUT_FILE}. "
                f"Found columns: {reader.fieldnames or []}",
                file=sys.stderr,
            )
            sys.exit(2)
        for i, row in enumerate(reader, start=2):   # row 1 is the header
            text = (row.get(required_column) or "").strip()
            comments.append((i, text))

    print(f"Analyzing {len(comments)} comments…\n")

    # ── Analyze ───────────────────────────────────────────────────────────────
    all_findings: List[Finding] = []
    rows_with_findings: List[Tuple[int, str, List[Finding]]] = []

    for row_num, text in comments:
        findings = analyze_comment(analyzer, text, row_num)
        rows_with_findings.append((row_num, text, findings))
        all_findings.extend(findings)
        if findings:
            summary = ", ".join(
                f"{f.entity_type}({f.score:.2f})" for f in findings
            )
            print(f"  Row {row_num:>4}: {len(findings)} finding(s) — {summary}")

    # ── Report ────────────────────────────────────────────────────────────────
    flagged = sum(1 for _, _, f in rows_with_findings if f)
    print(
        f"\nDone.  {len(all_findings)} finding(s) across "
        f"{flagged}/{len(comments)} comments.\n"
    )

    generate_html(rows_with_findings, all_findings, len(comments), OUTPUT_HTML)
    if all_findings:
        write_csv(all_findings, OUTPUT_CSV)
    else:
        print("  No findings — CSV not written.")


if __name__ == "__main__":
    main()

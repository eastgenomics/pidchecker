# PII Checker for Genomic Variant Classification Comments

A tool that scans free-text variant classification comments (exported as CSV) for personally identifiable information (PII) using [Microsoft Presidio](https://microsoft.github.io/presidio/), with custom recognizers tuned for the UK clinical genomics context.

## Why

Variant classification comments in genomic databases can inadvertently contain patient names, NHS numbers, dates of birth, postcodes, or internal patient/sample IDs. This tool flags those instances so they can be reviewed and redacted before data sharing or publication.

## What it detects

| Entity | Description |
|---|---|
| `PERSON` | Names (via spaCy NER) |
| `EMAIL_ADDRESS` | Email addresses |
| `URL` | URLs |
| `UK_NHS` | NHS numbers (with modulus-11 checksum validation) |
| `UK_POSTCODE` | UK postcodes |
| `DATE_OF_BIRTH` | Dates preceded by DOB context words |
| `PATIENT_ID` | Internal patient/sample IDs (e.g. `PAT-12345`, `SP12345R0123`, `GM1234567`) |

### False-positive suppression

Clinical genetics text is full of patterns that look like PII but aren't. The tool suppresses matches that overlap with:

- ACMG/AMP classification codes (e.g. `PM2`, `PVS1_Strong`)
- HGVS variant nomenclature (e.g. `c.1234A>G`, `NM_000527.4`)
- Literature citations (`Surname et al. 2021`)
- Genomic database names (gnomAD, ClinVar, OMIM, etc.)
- HGNC gene symbols (BRCA1, LDLR, CFTR, etc.)
- Journal abbreviations (J Med Genet, N Engl J Med, etc.)

## Setup

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python3 -m spacy download en_core_web_lg
```

## Usage

Place your input file as `comment.csv` in the project directory. It must contain a column named `comment_on_classification`.

```bash
.venv/bin/python3 pii_checker.py
```

### Output

- **`pii_report.html`** — interactive HTML report with highlighted PII spans (hover for entity type and confidence score)
- **`pii_findings.csv`** — machine-readable summary of all findings (row number, entity type, score, matched text, context snippet)

## Configuration

Edit the constants at the top of `pii_checker.py` to customise:

| Variable | Default | Description |
|---|---|---|
| `INPUT_FILE` | `comment.csv` | Path to the input CSV |
| `OUTPUT_HTML` | `pii_report.html` | Path for the HTML report |
| `OUTPUT_CSV` | `pii_findings.csv` | Path for the CSV findings |
| `SCORE_THRESHOLD` | `0.35` | Minimum confidence score (0.0–1.0). Lower = more findings, more false positives |
| `SPACY_MODEL` | `en_core_web_lg` | spaCy model (`en_core_web_lg` recommended; `en_core_web_trf` is more accurate but slower) |

# ClarVar: Open-Source Variant Annotation & Prioritisation Pipeline

---

## Overview

Clinical sequencing can produce **thousands of variants per sample**, but identifying the ones that actually matter is still largely a manual process — researchers chain together VEP, gnomAD lookups, and ClinVar queries by hand, introducing inconsistency and making workflows hard to reproduce.

ClarVar wraps that whole process into a single Python CLI tool. Give it a VCF file; it queries the Ensembl VEP REST API, pulls population allele frequencies from gnomAD, and cross-references ClinVar — then ranks your variants by a composite clinical priority score and writes the results to TSV, JSON, and a self-contained HTML report you can open in any browser.

It runs on a standard laptop. No local VEP installation, no database downloads, no HPC access required.

---

## Features

- Parses `.vcf` and `.vcf.gz` files; normalises chromosomes, skips structural variants
- Batched VEP REST API calls (up to 200 variants per request) with automatic rate-limit handling and exponential backoff
- Consequence annotation, CADD PHRED, SIFT, and PolyPhen-2 from a single API call
- Population allele frequencies from gnomAD v4.1 (730,947 exomes + 76,215 genomes)
- ClinVar clinical significance and variant IDs, somatic entries filtered out
- Evidence-based priority scoring (0–100) grounded in ACMG/AMP 2015 guidelines
- Three built-in ranking strategies: balanced default, rarity-first, clinical-evidence-first
- HPO gene filtering: restrict output to variants in genes associated with phenotype terms
- Four output formats: annotated VCF, ranked TSV, JSON summary, self-contained HTML report
- HTML report includes colour-coded ClinVar badges, priority score bars, and clickable OMIM and ClinVar links
- Clean Python API for embedding in custom pipelines
- 97 unit tests; opt-in live integration test against the Ensembl VEP API

---

## Quick Start

### Installation

```bash
git clone https://github.com/NosakhareOsaro/ClarVar.git
cd ClarVar
pip install -e ".[dev]"
```

Check everything is working:

```bash
clarvar --version
pytest tests/ -v        # 97 tests should pass
```

### Run the example

```bash
clarvar pipeline -i data/sample/sample_variants.vcf -o results/ --verbose
```

This writes four files to `results/`:

| File                       | Description                                                          |
| -------------------------- | -------------------------------------------------------------------- |
| `prioritised.vcf`          | Annotated VCF with CADD, gnomAD AF, and ClinVar in the INFO fields   |
| `prioritised_variants.tsv` | Ranked tab-separated table, one variant per row                      |
| `report.json`              | JSON summary with top variants and consequence/ClinVar distributions |
| `report.html`              | Self-contained HTML report — open in any browser, no internet needed |

### CLI reference

All commands write into the directory given by `-o`, which is created automatically if it doesn't exist.

```bash
# Full pipeline — annotation + prioritisation + all four output files
clarvar pipeline -i input.vcf -o results/

# Annotate only (writes annotated.vcf)
clarvar annotate -i raw.vcf -o results/

# Prioritise a pre-annotated VCF
clarvar prioritize -i annotated.vcf -o results/
```

**Available flags:**

| Flag              | Default  | Description                                                               |
| ----------------- | -------- | ------------------------------------------------------------------------- |
| `-i / --input`    | required | Input VCF file (`.vcf` or `.vcf.gz`)                                      |
| `-o / --output`   | required | Output directory (created if absent)                                      |
| `-a / --assembly` | `GRCh38` | Genome assembly — `GRCh38` or `GRCh37`                                    |
| `--no-html`       | off      | Skip HTML report generation                                               |
| `--top-n N`       | `50`     | Number of variants shown in the HTML report                               |
| `--hpo`           | none     | Filter output to genes associated with one or more HPO terms (repeatable) |
| `-v / --verbose`  | off      | Print step-by-step progress                                               |

```bash
# GRCh37 assembly
clarvar pipeline -i patient.vcf -o results/ -a GRCh37

# Skip HTML (useful in automated pipelines)
clarvar pipeline -i patient.vcf -o results/ --no-html

# Filter to variants in genes associated with seizures or intellectual disability
clarvar pipeline -i patient.vcf -o results/ --hpo HP:0001250 --hpo HP:0001249

# Show top 100 variants in the HTML report
clarvar pipeline -i patient.vcf -o results/ --top-n 100
```

---

## Python API

```python
from clarvar.vcf_parser import VCFParser
from clarvar.annotator import VariantAnnotator
from clarvar.prioritizer import VariantPrioritizer
from clarvar.reporter import write_tsv, write_html_report
from pathlib import Path

# 1. Parse
parser = VCFParser()
variants = parser.parse_file("patient.vcf")
print(f"Parsed {len(variants)} variants")

# 2. Annotate (calls Ensembl VEP REST API)
annotator = VariantAnnotator(assembly="GRCh38")
annotated = annotator.annotate_collection(variants, verbose=True)

# 3. Prioritise
prioritizer = VariantPrioritizer()
prioritized = prioritizer.prioritize_collection(annotated)

# 4. Output
write_tsv(prioritized, Path("results/variants.tsv"))
write_html_report(list(prioritized)[:50], Path("results/report.html"))

report = prioritizer.generate_report(prioritized)
print(f"Top variant: {report['top_variants'][0]}")
```

### Alternative ranking strategies

```python
from clarvar.prioritizer import RarityFirstStrategy, ClinicalEvidenceFirstStrategy

# Upweight population rarity — useful for rare disease workflows
ranked = RarityFirstStrategy().rank(list(annotated))

# Upweight ClinVar evidence — surfaces previously classified variants
ranked = ClinicalEvidenceFirstStrategy().rank(list(annotated))
```

### GRCh37 support

```python
annotator = VariantAnnotator(assembly="GRCh37")
```

### Offline / HPC annotation

For environments without internet access, use `LocalAnnotator` and supply your own pre-downloaded gnomAD and ClinVar databases:

```python
from clarvar.annotator import LocalAnnotator

la = LocalAnnotator(db_path="/data/annotation_db")
la.load_databases()
annotated = la.annotate_collection(variants)
```

### HPO gene filtering

```python
from clarvar.hpo_filter import filter_collection_by_hpo

# Restrict to variants in genes associated with the given HPO terms
# The annotation file is downloaded and cached automatically on first use
filtered = filter_collection_by_hpo(prioritized, ["HP:0001250", "HP:0001249"])
```

---

## How It Works

### Pipeline

```
Input VCF (.vcf / .vcf.gz)
        │
        ▼
[VCFParser]  ── parse variants, normalise chromosomes, skip SVs
        │
        ▼
[VariantAnnotator]  ── batch POST to Ensembl VEP REST API
        │                 ├── consequence + HGVS notation
        │                 ├── CADD PHRED / SIFT / PolyPhen-2
        │                 ├── gnomAD v4.1 allele frequency
        │                 └── ClinVar clinical significance
        ▼
[VariantPrioritizer]  ── score + sort
        │
        ▼
[Reporter]  ── TSV table  +  HTML report  +  JSON summary  +  annotated VCF
```

### Priority score (0–100)

The scoring model is grounded in the ACMG/AMP 2015 variant classification guidelines. Consequence and clinical evidence are weighted equally as primary criteria; population frequency is a supporting criterion.

| Component            | Weight | Example values                                                        |
| -------------------- | ------ | --------------------------------------------------------------------- |
| Consequence severity | 40%    | Frameshift/stop = 100, missense = 40, synonymous = 10, intergenic = 0 |
| ClinVar significance | 40%    | Pathogenic = 100, likely pathogenic = 80, VUS = 30, benign = 5        |
| Population rarity    | 20%    | AF < 0.01% → 100, AF 0.01–1% → 80, AF > 10% → 20                      |

**Final score** = 0.4 × consequence + 0.4 × ClinVar + 0.2 × rarity

When no gnomAD frequency is available, the rarity component defaults to 50 — absence of population data is not the same as being common. When no ClinVar entry exists, the ClinVar component defaults to 20.

---

## Data Sources

| Source                                           | What it provides                        | Version                                |
| ------------------------------------------------ | --------------------------------------- | -------------------------------------- |
| [Ensembl VEP REST](https://rest.ensembl.org)     | Consequence, HGVS, CADD, SIFT, PolyPhen | v115 (GRCh38), v111 (GRCh37)           |
| [gnomAD](https://gnomad.broadinstitute.org)      | Population allele frequencies           | v4.1 (730,947 exomes, 76,215 genomes)  |
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Clinical significance classifications   | Weekly release                         |
| [HPO](https://hpo.jax.org)                       | Gene-to-phenotype associations          | genes_to_phenotype.txt, cached locally |

All sources are publicly available and require no data access agreement.

**References:**

1. McLaren W et al. (2016) The Ensembl Variant Effect Predictor. _Genome Biology_, 17:122.
2. Karczewski KJ et al. (2020) The mutational constraint spectrum from variation in 141,456 humans. _Nature_, 581:434–443.
3. Landrum MJ et al. (2018) ClinVar: improving access to variant interpretations. _Nucleic Acids Research_, 46:D1062–D1067.
4. Rentzsch P et al. (2019) CADD: predicting the deleteriousness of variants throughout the human genome. _Nucleic Acids Research_, 47:D886–D894.
5. Richards S et al. (2015) Standards and guidelines for the interpretation of sequence variants. _Genetics in Medicine_, 17:405–423.

---

## Project Structure

```
ClarVar/
├── src/clarvar/
│   ├── __init__.py
│   ├── variant.py        # Variant + VariantCollection dataclasses
│   ├── vcf_parser.py     # VCF I/O (.vcf and .vcf.gz)
│   ├── annotator.py      # Ensembl VEP REST API integration
│   ├── prioritizer.py    # Scoring, filtering, ranking strategies
│   ├── reporter.py       # TSV + HTML output
│   ├── hpo_filter.py     # HPO gene filtering
│   └── cli.py            # click-based CLI
├── workflow/
│   ├── Snakefile         # Snakemake workflow for cohort-scale runs
│   ├── config.yaml       # Input/output paths and assembly config
│   └── README.md         # Workflow documentation
├── tests/
│   ├── conftest.py       # Shared fixtures
│   ├── test_variant.py
│   ├── test_annotator.py # Mocked API tests + optional live tests
│   ├── test_prioritizer.py
│   └── test_reporter.py
├── data/sample/
│   └── sample_variants.vcf
├── pyproject.toml
├── CONTRIBUTING.md
├── LICENSE
└── README.md
```

---

## Running Tests

```bash
# Full offline test suite (97 tests, no API calls)
pytest tests/ -v

# Include live Ensembl VEP API tests (requires internet)
CLARVAR_LIVE_TESTS=1 pytest tests/ -v

# With coverage report
pytest tests/ --cov=clarvar --cov-report=html
```

---

## Snakemake Workflow

For cohort-scale annotation across multiple samples, a Snakemake workflow is provided in `workflow/`. It supports both directory-based VCF discovery and explicit sample sheets.

```bash
# Run on all VCFs in the configured input directory
snakemake -s workflow/Snakefile --cores 4

# With a sample sheet (sample,vcf CSV format)
# Set sample_sheet in workflow/config.yaml, then:
snakemake -s workflow/Snakefile --cores 8
```

Per-sample outputs are written to `results/<sample>/`. A cohort-level summary is written to `results/cohort_summary.json`.

---

## Requirements

- Python ≥ 3.10
- `click` ≥ 8.1
- `requests` ≥ 2.31
- Internet access for VEP annotation (or use `LocalAnnotator` for offline mode)
- Snakemake (optional, for the cohort workflow only)

---

## Comparison with Other Tools

| Feature                    | ClarVar | bcftools | Ensembl VEP |
| -------------------------- | ------- | -------- | ----------- |
| Lightweight (pip install)  | ✅      | ✅       | ❌          |
| Real VEP annotation        | ✅      | ❌       | ✅          |
| gnomAD v4.1 frequencies    | ✅      | ✅       | ✅          |
| ClinVar integration        | ✅      | ✅       | ✅          |
| CADD / SIFT / PolyPhen     | ✅      | ❌       | ✅          |
| Composite priority scoring | ✅      | ❌       | ❌          |
| HTML report with links     | ✅      | ❌       | ❌          |
| HPO gene filtering         | ✅      | ❌       | ❌          |
| Snakemake cohort workflow  | ✅      | ❌       | ❌          |
| Python API                 | ✅      | ❌       | ✅          |
| No local DB required       | ✅      | ❌       | ❌          |

---

## Limitations

- Requires internet access for VEP annotation — use `LocalAnnotator` for offline/HPC use
- Annotation rate is approximately 200 variants per 3–5 seconds (Ensembl VEP rate limit); large cohorts benefit from local VEP or the Snakemake parallel workflow
- Single-sample VCF files only; multi-sample support is planned
- Structural variants are currently skipped
- Error handling and input validation are minimal in v0.1.0 — this is the first priority for v0.2.0
- The priority score is a research tool for variant triage, not a clinical diagnostic

---

## Roadmap

**v0.1.0** (current) — Full annotation pipeline, 97 tests, HPO gene filtering, Snakemake cohort workflow, four output formats

**v0.2.0** — Comprehensive error handling and input validation throughout; better error messages for common failure modes

**v0.3.0** — Local VEP offline mode; Docker/Singularity wrapper for HPC environments

**v0.4.0** — Multi-sample VCF support; ACMG criteria scoring

**v1.0.0** — Benchmarked against ClinVar pathogenic/benign gold-standard set with precision/recall evaluation; performance-optimised for large cohorts

---

## Contributing

Contributions are welcome. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for the full workflow, coding standards, and PR guidelines.

```bash
# Fork and clone
git clone https://github.com/YOUR-USERNAME/ClarVar.git
cd ClarVar

# Create a branch
git checkout -b feature/my-feature

# Install in dev mode
pip install -e ".[dev]"

# Make changes, write tests, then
pytest tests/ -v
git commit -m "feat: description of change"
git push origin feature/my-feature
# → open a Pull Request on GitHub
```

Issues are labelled by difficulty — look for **good first issue** to get started.

---

## Citation

If you use ClarVar in your research, please cite:

```bibtex
@software{clarvar2026,
  author       = {Osaro, Nosakhare Odionfo and
                  Rodriguez Sosa, Alejandra and
                  Gimenez-Rios, Amparo and
                  Darling, Benjamin and
                  Biju Daniel, Abel and
                  Ravichandran, Lochan},
  title        = {ClarVar: Open-Source Variant Annotation and Prioritisation Pipeline},
  year         = {2026},
  version      = {0.1.0},
  url          = {https://github.com/NosakhareOsaro/ClarVar},
  license      = {MIT}
}
```

---

## License

MIT — free to use, modify, and distribute. See [LICENSE](LICENSE) for details.

---

## Acknowledgements

Built during the University of Glasgow CompBio Hackathon 2026 by a team of MSc Bioinformatics students and researchers. Data from Ensembl VEP, gnomAD, ClinVar, and the Human Phenotype Ontology. Inspired by the nf-core community's commitment to reproducible, well-documented bioinformatics tooling.

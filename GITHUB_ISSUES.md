# ClarVar — GitHub Issues

Paste each block below directly into a new GitHub Issue.
Go to: https://github.com/NosakhareOsaro/ClarVar/issues/new

---

## Issue 1 · good first issue

**Title:** Convert a Variant to Ensembl VEP region notation

**Labels:** `good first issue`, `annotator`

**Body:**

The Ensembl VEP REST API requires variants to be submitted in a specific region
string format. Implement the `_to_vep_region()` helper in `src/clarvar/annotator.py`.

**Format:** `"{chrom} {pos} {end} {ref}/{alt} 1"`
where `end = pos + len(ref) - 1`

**Examples:**
- SNV `1:100 A>T` → `"1 100 100 A/T 1"`
- Deletion `7:100 ATG>A` → `"7 100 102 ATG/A 1"`
- Insertion `2:100 A>ATG` → `"2 100 100 A/ATG 1"`

**Acceptance criteria:**
- All three `TestVepRegionNotation` tests in `tests/test_annotator.py` pass.

**File:** `src/clarvar/annotator.py` → `_to_vep_region()`

---

## Issue 2 · good first issue

**Title:** Extract gnomAD allele frequency from a VEP response

**Labels:** `good first issue`, `annotator`

**Body:**

Implement `_extract_gnomad_af()` in `src/clarvar/annotator.py`.

gnomAD frequencies are nested inside the VEP response under:
`hit["colocated_variants"][i]["frequencies"][allele]["gnomad"]`

Also check for keys `"gnomadg"` (genomes) and `"gnomade"` (exomes) as fallbacks.
Return the first value found, or `None` if absent.

**Acceptance criteria:**
- All five `TestGnomadExtraction` tests in `tests/test_annotator.py` pass.

**File:** `src/clarvar/annotator.py` → `_extract_gnomad_af()`

---

## Issue 3 · good first issue

**Title:** Extract ClinVar significance from a VEP response

**Labels:** `good first issue`, `annotator`

**Body:**

Implement `_extract_clinvar()` in `src/clarvar/annotator.py`.

ClinVar data lives in `hit["colocated_variants"][i]["clin_sig"]` (a list).
The ClinVar variation ID lives in `hit["colocated_variants"][i]["var_synonyms"]["ClinVar"]`.

Rules:
- Skip entries where `"somatic": True`.
- Join multiple significance values with `", "`.
- Return `(None, None)` if no ClinVar data is found.

**Acceptance criteria:**
- All four `TestClinvarExtraction` tests in `tests/test_annotator.py` pass.

**File:** `src/clarvar/annotator.py` → `_extract_clinvar()`

---

## Issue 4 · good first issue

**Title:** Select the most severe VEP transcript consequence

**Labels:** `good first issue`, `annotator`

**Body:**

Implement `_pick_most_severe_transcript()` in `src/clarvar/annotator.py`.

VEP returns a list of `transcript_consequences`, each with a `"consequence_terms"` list.
Use the `CONSEQUENCE_SEVERITY` dict (already defined) to rank terms — lower rank = more severe.
Return the full transcript dict (not just the term) so the caller can access `gene_symbol`,
`hgvsc`, `cadd_phred`, etc.

Return `None` if the list is empty.

**File:** `src/clarvar/annotator.py` → `_pick_most_severe_transcript()`

---

## Issue 5 · intermediate

**Title:** Apply a VEP result to a Variant object

**Labels:** `intermediate`, `annotator`

**Body:**

Implement `_apply_vep_result()` in `src/clarvar/annotator.py`.

This function ties together the three helper functions above and writes
all annotation fields onto a `Variant` object in-place.

Fields to populate:
- From the most severe transcript: `gene`, `transcript`, `consequence`,
  `hgvsc`, `hgvsp`, `sift`, `polyphen`, `cadd_phred`
- From colocated_variants: `allele_frequency`, `clinvar_significance`, `clinvar_id`

**Acceptance criteria:**
- All seven `TestApplyVepResult` tests in `tests/test_annotator.py` pass.

**File:** `src/clarvar/annotator.py` → `_apply_vep_result()`

---

## Issue 6 · intermediate

**Title:** POST a batch of variants to the Ensembl VEP REST API

**Labels:** `intermediate`, `annotator`

**Body:**

Implement `_post_vep_batch()` in `src/clarvar/annotator.py`.

The function should:
1. Build a JSON payload with `"variants"`, `"CADD": 1`, `"gnomAD": 1`, `"ClinVar": 1`,
   `"canonical": 1`, `"pick": 1`.
2. POST to the given endpoint with `Content-Type: application/json`.
3. On HTTP 200, return the parsed JSON list.
4. On HTTP 429 (rate limited), wait `Retry-After` seconds and retry.
5. On other errors, retry up to `MAX_RETRIES` times with exponential backoff,
   then log a warning and return `[]`.

**Files:** `src/clarvar/annotator.py` → `_post_vep_batch()`

**Testing tip:** Use `unittest.mock.patch` to mock `requests.post` — no real
API calls needed in unit tests. See `TestVariantAnnotatorMocked` in `test_annotator.py`.

---

## Issue 7 · intermediate

**Title:** Implement VariantAnnotator.annotate_collection()

**Labels:** `intermediate`, `annotator`

**Body:**

Implement `VariantAnnotator.annotate_collection()` in `src/clarvar/annotator.py`.

This is the main public method. It should:
1. Build a lookup dict: `{vep_region_string: Variant}` using `_to_vep_region()`.
2. Send variants in batches of `BATCH_SIZE` using `_post_vep_batch()`.
3. For each hit in the response, look up the matching Variant via `hit["input"]`
   and call `_apply_vep_result()`.
4. Sleep `RATE_LIMIT_PAUSE` seconds between batches.
5. Return a new `VariantCollection` with the annotated variants.

**Acceptance criteria:**
- `TestVariantAnnotatorMocked` tests in `tests/test_annotator.py` pass.

**Files:** `src/clarvar/annotator.py` → `VariantAnnotator.annotate_collection()`

---

## Issue 8 · good first issue

**Title:** Implement variant priority scoring

**Labels:** `good first issue`, `prioritizer`

**Body:**

Implement the three scoring methods in `src/clarvar/prioritizer.py`:
- `_score_consequence()` — uses `CONSEQUENCE_SCORES` dict
- `_score_clinvar()` — uses `CLINVAR_SCORES` dict
- `_score_frequency()` — uses the frequency tier logic in the module docstring

Also populate the `CONSEQUENCE_SCORES` and `CLINVAR_SCORES` class dicts
with the values from the module docstring.

**Acceptance criteria:**
- `TestVariantPrioritizer::test_score_consequence`,
  `test_score_clinvar`, `test_score_frequency` all pass.

**File:** `src/clarvar/prioritizer.py`

---

## Issue 9 · good first issue

**Title:** Implement VariantPrioritizer.prioritize_collection()

**Labels:** `good first issue`, `prioritizer`

**Body:**

Implement `prioritize_variant()` and `prioritize_collection()` in
`src/clarvar/prioritizer.py`.

`prioritize_variant()` computes:
```
score = (consequence_weight × _score_consequence())
      + (clinvar_weight     × _score_clinvar())
      + (frequency_weight   × _score_frequency())
```
and stores it in `variant.priority_score`.

`prioritize_collection()` calls `prioritize_variant()` on each variant,
then returns a new `VariantCollection` sorted descending by `priority_score`.

**Acceptance criteria:**
- `TestVariantPrioritizer::test_prioritize_variant` and
  `test_prioritize_collection` pass.

**File:** `src/clarvar/prioritizer.py`

---

## Issue 10 · intermediate

**Title:** Implement VariantPrioritizer.generate_report() and ranking strategies

**Labels:** `intermediate`, `prioritizer`

**Body:**

Implement `generate_report()`, `_get_consequence_distribution()`, and
`_get_clinvar_distribution()` in `src/clarvar/prioritizer.py`.

Also implement the three `RankingStrategy` subclasses:
- `ConsequenceOnlyStrategy` — consequence weight 1.0, others 0.0
- `RarityFirstStrategy` — frequency weight 0.6, others 0.2
- `ClinicalEvidenceFirstStrategy` — clinvar weight 0.6, others 0.2

**Acceptance criteria:**
- `TestVariantPrioritizer::test_generate_report`,
  `test_consequence_distribution`, `test_clinvar_distribution` pass.
- All `TestRankingStrategies` tests pass.

**File:** `src/clarvar/prioritizer.py`

---

## Issue 11 · good first issue

**Title:** Implement write_tsv() — export variants to a TSV file

**Labels:** `good first issue`, `reporter`

**Body:**

Implement `write_tsv()` and `_fmt()` in `src/clarvar/reporter.py`.

`_fmt()` should:
- Return `"."` for `None`
- Format floats to `precision` decimal places
- Return `str(value)` for everything else

`write_tsv()` should write a tab-separated file with a header row matching
`TSV_COLUMNS` and one data row per variant, with `None` values as `"."`.

**Acceptance criteria:**
- All `TestWriteTsv` and `TestFmt` tests in `tests/test_reporter.py` pass.

**File:** `src/clarvar/reporter.py` → `write_tsv()`, `_fmt()`

---

## Issue 12 · intermediate

**Title:** Implement write_html_report() — self-contained HTML report

**Labels:** `intermediate`, `reporter`

**Body:**

Implement `write_html_report()`, `_clinvar_badge()`, and `_priority_bar()`
in `src/clarvar/reporter.py`.

The HTML report must:
- Be a single self-contained file (no CDN links — must work offline)
- Include a header with filename and timestamp
- Show a filter summary bar
- Render a table with colour-coded ClinVar badges and score bars
- Link gene names to OMIM and ClinVar IDs to ncbi.nlm.nih.gov/clinvar

All user-supplied strings must go through `html.escape()` before insertion.

**Acceptance criteria:**
- All `TestWriteHtmlReport` tests in `tests/test_reporter.py` pass.

**File:** `src/clarvar/reporter.py`

---

## Issue 13 · intermediate

**Title:** Implement the `clarvar annotate` CLI command

**Labels:** `intermediate`, `cli`

**Body:**

Implement the `annotate` click command in `src/clarvar/cli.py`.

The command should:
1. Accept `-i` (input VCF), `-o` (output dir), `-a` (assembly), `-v` (verbose)
2. Create the output directory
3. Parse the VCF, annotate with `VariantAnnotator`, write `<output>/annotated.vcf`
4. Print progress with `click.echo` and success with `click.secho(..., fg="green")`

**File:** `src/clarvar/cli.py` → `annotate()`

---

## Issue 14 · intermediate

**Title:** Implement the `clarvar pipeline` CLI command

**Labels:** `intermediate`, `cli`

**Body:**

Implement the `pipeline` click command in `src/clarvar/cli.py` — the main
entry point that runs annotation + prioritisation in one step.

Output files written to the `-o` directory:
- `prioritised.vcf`
- `prioritised_variants.tsv`
- `report.json`
- `report.html` (unless `--no-html`)

Additional flags: `--assembly`, `--no-html`, `--top-n N`

**File:** `src/clarvar/cli.py` → `pipeline()`

---

## Issue 15 · stretch

**Title:** HPO gene filtering — restrict output to phenotype-relevant genes

**Labels:** `stretch`, `annotator`

**Body:**

Add an optional `--hpo` flag to the `clarvar pipeline` command that accepts
a list of Human Phenotype Ontology (HPO) term IDs (e.g. `HP:0001250`).

When provided, only variants in genes associated with those HPO terms should
be included in the prioritised output.

**Data source:** The HPO gene-to-phenotype annotations are freely available at
https://hpo.jax.org/data/annotations

**Suggested approach:**
1. Download/cache the `genes_to_phenotype.txt` file
2. Build a lookup: HPO term → set of gene symbols
3. Filter the prioritised VariantCollection before writing outputs

**File:** New `src/clarvar/hpo_filter.py` + integration in `cli.py`

---

## Issue 16 · stretch

**Title:** Snakemake workflow for cohort-scale annotation

**Labels:** `stretch`, `workflow`

**Body:**

Create a `Snakefile` that wraps the `clarvar pipeline` command to support
parallel annotation of multiple VCF files.

The workflow should:
1. Accept a directory of input VCF files (or a sample sheet)
2. Run `clarvar pipeline` for each sample in parallel
3. Produce per-sample output directories
4. Generate a combined summary report across all samples

**File:** `workflow/Snakefile` + `workflow/config.yaml`

**Reference:** https://snakemake.readthedocs.io

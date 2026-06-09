# ClarVar Architecture

This document describes the design and architecture of ClarVar.

## Table of Contents

- [Overview](#overview)
- [System Design](#system-design)
- [Data Flow](#data-flow)
- [Module Overview](#module-overview)
- [Key Concepts](#key-concepts)
- [Extensibility](#extensibility)

---

## Overview

ClarVar follows a **pipeline architecture** with three main stages:

```
Input VCF → [PARSE] → [ANNOTATE] → [PRIORITIZE] → Output VCF + Report
```

Each stage can be run independently or as part of the complete pipeline.

### Design Principles

1. **Modularity** - Each component can be used independently
2. **Simplicity** - Minimal dependencies, easy to understand
3. **Reproducibility** - Deterministic, version-controlled results
4. **Extensibility** - Easy to add new annotation sources or ranking strategies
5. **Robustness** - Graceful error handling, informative logging

---

## System Design

### High-Level Architecture

```
┌─────────────────────────────────────────────────────┐
│                  CLI Interface                       │
│  (clarvar annotate/prioritize/pipeline commands)    │
└────────────────┬────────────────────────────────────┘
                 │
        ┌────────┴─────────┬──────────────┐
        ▼                  ▼              ▼
    ┌────────┐         ┌──────────┐  ┌──────────┐
    │  VCF   │         │ Python   │  │   JSON   │
    │ Parser │         │   API    │  │  Output  │
    └────────┘         └──────────┘  └──────────┘
        │                  │              │
        └──────────┬───────┴──────────────┘
                   ▼
        ┌──────────────────────┐
        │  Core Processing     │
        │  Engine              │
        ├──────────────────────┤
        │ VariantCollection    │
        │ VariantAnnotator     │
        │ VariantPrioritizer   │
        └──────────────────────┘
                   │
        ┌──────────┴──────────────┐
        ▼                         ▼
    ┌──────────┐           ┌─────────────┐
    │ External │           │   Local     │
    │   APIs   │           │ Databases   │
    └──────────┘           └─────────────┘
```

### Module Organization

```
clarvar/
├── variant.py          # Variant data model
├── annotator.py        # Annotation engines
├── prioritizer.py      # Prioritization algorithms
├── vcf_parser.py       # File I/O (VCF format)
├── cli.py              # Command-line interface
└── __init__.py         # Package exports
```

---

## Data Flow

### Pipeline Stage 1: Parsing

```
Raw VCF File
    │
    ├─ Read file line by line
    ├─ Skip comment lines (#)
    ├─ Parse CHROM POS REF ALT columns
    ├─ Extract optional fields (QUAL, FORMAT, GENOTYPE)
    │
    ▼
Variant Objects
    │
    ├─ Create Variant dataclass instances
    ├─ Validate chromosome and position
    ├─ Store in VariantCollection
    │
    ▼
VariantCollection (in memory)
```

### Pipeline Stage 2: Annotation

```
VariantCollection
    │
    ├─ For each Variant:
    │  ├─ Query gnomAD for population frequency
    │  ├─ Query ClinVar for clinical significance
    │  ├─ Predict functional consequence
    │  └─ Add annotations to Variant object
    │
    ▼
Annotated VariantCollection
    │
    ├─ Each variant now has:
    │  ├─ allele_frequency
    │  ├─ clinvar_significance
    │  └─ consequence
    │
    ▼
Ready for prioritization
```

### Pipeline Stage 3: Prioritization

```
Annotated VariantCollection
    │
    ├─ For each Variant:
    │  ├─ Score consequence (0-100)
    │  ├─ Score ClinVar significance (0-100)
    │  ├─ Score population frequency (0-100)
    │  ├─ Combine with weighted formula:
    │  │  priority = 0.4×consequence + 0.4×clinvar + 0.2×frequency
    │  └─ Store priority_score in variant
    │
    ├─ Sort by priority_score (descending)
    │
    ▼
Prioritized VariantCollection
    │
    ├─ Save to VCF with annotations
    ├─ Generate JSON report
    │
    ▼
Results (VCF + JSON)
```

---

## Module Overview

### 1. Variant Module (`variant.py`)

**Purpose:** Define variant data structures

**Key Classes:**

- **`Variant`** - Single genetic variant
  - Attributes: chromosome, position, ref, alt, annotations
  - Methods: type inference, rarity checking, pathogenicity assessment
  
- **`VariantCollection`** - Container for multiple variants
  - Attributes: variants list, collection name
  - Methods: filter, sort, count, top_n, iteration

**Key Enums:**

- `VariantType` - SNP, INDEL, DELETION, INSERTION, STRUCTURAL
- `Consequence` - FRAMESHIFT, STOP_GAINED, MISSENSE, SYNONYMOUS, etc.

**Design Notes:**

- Variants are immutable after creation (good for data integrity)
- Collections support filtering and sorting operations
- All operations return new collections (functional style)

### 2. Annotator Module (`annotator.py`)

**Purpose:** Enrich variants with external data

**Key Classes:**

- **`VariantAnnotator`** - Queries external APIs
  - Methods:
    - `annotate_variant()` - Annotate single variant
    - `annotate_collection()` - Annotate all variants
    - `_add_gnomad_frequency()` - Add population frequency
    - `_add_clinvar_info()` - Add clinical significance
    - `_predict_consequence()` - Predict functional impact

- **`LocalAnnotator`** - Uses pre-downloaded databases
  - Methods:
    - `load_databases()` - Load from local files
    - `annotate_variant()` - Annotate from cache

**API Integration:**

- gnomAD API: Population allele frequencies
- ClinVar API: Clinical variant classifications
- Ensembl REST API: Consequence prediction (placeholder)

**Design Notes:**

- Currently uses simulated data for development
- Production version queries real APIs with caching
- Graceful error handling if APIs are unavailable
- Supports both online and offline modes

### 3. Prioritizer Module (`prioritizer.py`)

**Purpose:** Rank variants by clinical importance

**Key Classes:**

- **`VariantPrioritizer`** - Weighted scoring system
  - Attributes: consequence_weight, clinvar_weight, frequency_weight
  - Methods:
    - `prioritize_variant()` - Score single variant
    - `prioritize_collection()` - Score and sort all variants
    - `generate_report()` - Create summary statistics
    - `_score_consequence()`, `_score_clinvar()`, `_score_frequency()`

- **`RankingStrategy`** - Base class for different strategies
  - `ConsequenceOnlyStrategy` - Rank by impact alone
  - `RarityFirstStrategy` - Prioritize rare variants
  - `ClinicalEvidenceFirstStrategy` - Prioritize known pathogenic variants

**Scoring Algorithm:**

```
priority_score = 0.4 × consequence_score +
                 0.4 × clinvar_score +
                 0.2 × frequency_score

Where:
  consequence_score ∈ [0, 100]  (frameshift=100, missense=40, synonymous=10)
  clinvar_score ∈ [0, 100]       (pathogenic=100, benign=5)
  frequency_score ∈ [0, 100]     (rare=100, common=20)
```

**Design Notes:**

- Weights are normalized so they sum to 1.0
- Modular scoring functions allow easy customization
- Multiple ranking strategies support different use cases
- Reports include distribution statistics and top variants

### 4. VCF Parser Module (`vcf_parser.py`)

**Purpose:** Read and write VCF format files

**Key Classes:**

- **`VCFParser`** - VCF format handler
  - Methods:
    - `parse_file()` - Read VCF file
    - `parse_lines()` - Parse VCF format lines
    - `write_vcf()` - Write annotated VCF
    - `_parse_variant_line()` - Parse single variant line
    - `_parse_header_line()` - Extract sample information

**VCF Format Support:**

- Reads VCFv4.2 standard format
- Handles variable columns and optional fields
- Writes enriched VCF with annotation data in INFO field
- Preserves sample genotype information

**Helper Functions:**

- `parse_simple_tsv()` - Parse tab-separated variant files

**Design Notes:**

- Lazy parsing (one variant at a time)
- Memory efficient for large files
- Comprehensive error handling
- Logging for debugging

### 5. CLI Module (`cli.py`)

**Purpose:** Command-line interface

**Key Classes:**

- **`CLIInterface`** - CLI command handler
  - Methods:
    - `main()` - Entry point
    - `annotate_command()` - CLI: annotation only
    - `prioritize_command()` - CLI: prioritization only
    - `pipeline_command()` - CLI: full pipeline

**Commands:**

```
clarvar annotate    - Annotate variants
clarvar prioritize  - Prioritize variants
clarvar pipeline    - Run complete pipeline
```

**Design Notes:**

- Works without click library (fallback to argparse)
- Comprehensive logging and progress reporting
- User-friendly error messages
- Supports both interactive and scripted usage

---

## Key Concepts

### Variant Representation

A variant is defined by its location and alleles:

```
Chromosome:Position Reference>Alternate

Example: 17:41246881 G>A
```

### Consequence Hierarchy (by severity)

```
100  ├─ Frameshift
     ├─ Stop gained
     ├─ Stop lost
     ├─ Splice acceptor/donor
 80  │
 60  ├─ Missense
 40  │
 20  ├─ Upstream/Downstream
 10  ├─ Synonymous
  0  └─ Intergenic
```

### Allele Frequency Classes

```
< 0.0001    Very rare (new mutations)
0.0001-0.01 Rare (likely pathogenic if functional)
0.01-0.05   Low frequency
0.05-0.5    Common
> 0.5       Very common (likely benign)
```

### Filtering Chain

```
All Variants
    ↓
[Filter by consequence] → High-impact variants
    ↓
[Filter by rarity] → Rare high-impact variants
    ↓
[Filter by ClinVar] → Known disease variants
    ↓
[Sort by score] → Prioritized list
```

---

## Data Structures

### Variant Object

```python
@dataclass
class Variant:
    # Identification
    chromosome: str
    position: int
    ref: str              # Reference allele
    alt: str              # Alternate allele
    sample_id: str        # Optional
    
    # Quality
    quality: float        # VCF QUAL field
    genotype: str         # e.g., "0/1" (heterozygous)
    
    # Annotation Results
    consequence: Consequence         # Functional impact
    gene: Optional[str]              # Associated gene
    allele_frequency: Optional[float]
    clinvar_significance: Optional[str]
    clinvar_conditions: List[str]
    
    # Scoring
    priority_score: float            # Final score [0-100]
```

### VariantCollection Object

```python
@dataclass
class VariantCollection:
    variants: List[Variant]          # Collection of variants
    name: str                        # Collection name
    
    # Methods
    add(variant)                     # Add single variant
    filter_by_consequence()          # Filter by impact
    filter_by_rarity()              # Filter by frequency
    sort_by_priority()              # Sort by score
    top_n(n)                        # Get top N variants
    count_by_gene()                 # Gene statistics
```

---

## Extensibility

### Adding a New Annotation Source

```python
class MyAnnotator:
    def annotate_variant(self, variant: Variant) -> Variant:
        # 1. Query your data source
        my_data = self.query_my_api(variant)
        
        # 2. Update variant object
        variant.annotation['my_field'] = my_data
        
        # 3. Return updated variant
        return variant
```

### Adding a New Ranking Strategy

```python
from clarvar.prioritizer import RankingStrategy

class MyStrategy(RankingStrategy):
    def rank(self, variants: List[Variant]) -> List[Variant]:
        # Custom ranking logic
        ranked = sorted(
            variants,
            key=lambda v: my_custom_score(v),
            reverse=True
        )
        return ranked
```

### Custom Filtering

```python
# Filter variants with custom logic
filtered = [
    v for v in collection
    if v.priority_score > 50 and v.allele_frequency < 0.01
]
```

---

## Performance Characteristics

### Time Complexity

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| Parse N variants | O(N) | Single pass through file |
| Annotate N variants | O(N) | Serial processing |
| Prioritize N variants | O(N) | Single pass scoring |
| Filter N variants | O(N) | Linear filtering |
| Sort N variants | O(N log N) | Standard sort |

### Space Complexity

| Data Structure | Complexity | Notes |
|----------------|-----------|-------|
| VariantCollection | O(N) | One Variant per record |
| Single Variant | O(1) | Fixed-size + annotations |
| Annotation Cache | O(M) | M = number of cached API calls |

### Performance Tips

1. **Process in batches** - Split large files
2. **Filter early** - Reduce data size before expensive operations
3. **Use local databases** - Avoid repeated API calls
4. **Enable caching** - Reuse annotation results

---

## Error Handling

### Graceful Degradation

```
┌─ API available
├─ API timeout → Use defaults
├─ Invalid file → Skip record
├─ Missing data → Continue with warnings
└─ Critical error → Log and fail gracefully
```

### Logging Levels

- **DEBUG**: Detailed processing information
- **INFO**: High-level progress updates
- **WARNING**: Issues that don't stop processing
- **ERROR**: Failures that stop a step
- **CRITICAL**: System failures

---

## Testing Strategy

### Unit Tests

- Variant creation and validation
- Collection operations (filter, sort, etc.)
- Annotation logic
- Prioritization scoring

### Integration Tests

- End-to-end pipeline
- File I/O (parsing, writing)
- CLI commands

### Test Coverage Goals

- 80%+ code coverage
- All public APIs tested
- Error paths tested

---

## Future Enhancements

### Short Term (v0.2)

- Real API integration
- Local database support
- Performance optimization
- More ranking strategies

### Medium Term (v0.3)

- Multi-sample VCF support
- Web interface
- Variant visualization
- Batch processing tools

### Long Term (v1.0)

- Structural variant support
- Machine learning ranking
- Clinical phenotype integration
- Integrating with EHR systems

---

## References

- VCF Specification: https://samtools.github.io/hts-specs/VCFv4.2.pdf
- gnomAD: https://gnomad.broadinstitute.org/
- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/
- Ensembl: https://www.ensembl.org/

---

**Last Updated:** June 2026  
**Version:** 0.1.0

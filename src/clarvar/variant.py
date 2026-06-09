"""
Variant class for representing and managing genetic variants.
"""

from dataclasses import dataclass, field, asdict
from typing import Optional, List, Dict, Any
from enum import Enum


class VariantType(Enum):
    """Enumeration of variant types."""
    SNP = "SNP"
    INDEL = "INDEL"
    DELETION = "DELETION"
    INSERTION = "INSERTION"
    STRUCTURAL = "STRUCTURAL"
    UNKNOWN = "UNKNOWN"


class Consequence(Enum):
    """
    Enumeration of VEP consequence terms ordered by severity.

    Note: only one value per consequence string — the duplicate
    MISSENSE_VARIANT entry that caused silent scoring bugs has been removed.
    """
    TRANSCRIPT_ABLATION = "transcript_ablation"
    SPLICE_ACCEPTOR = "splice_acceptor_variant"
    SPLICE_DONOR = "splice_donor_variant"
    STOP_GAINED = "stop_gained"
    FRAMESHIFT = "frameshift_variant"
    STOP_LOST = "stop_lost"
    START_LOST = "start_lost"
    TRANSCRIPT_AMPLIFICATION = "transcript_amplification"
    INFRAME_INSERTION = "inframe_insertion"
    INFRAME_DELETION = "inframe_deletion"
    DISRUPTIVE_INFRAME = "disruptive_inframe_deletion"
    MISSENSE = "missense_variant"
    SPLICE_REGION = "splice_region_variant"
    SYNONYMOUS = "synonymous_variant"
    UPSTREAM = "upstream_gene_variant"
    DOWNSTREAM = "downstream_gene_variant"
    INTRON = "intron_variant"
    INTERGENIC = "intergenic_variant"
    UNKNOWN = "unknown"

    @classmethod
    def from_vep_string(cls, term: str) -> "Consequence":
        """Look up a Consequence by its VEP term string, defaulting to UNKNOWN."""
        for member in cls:
            if member.value == term:
                return member
        return cls.UNKNOWN


@dataclass
class Variant:
    """
    Represents a single genetic variant.

    Attributes
    ----------
    chromosome : str
        Chromosome (1-22, X, Y, MT). No 'chr' prefix — normalised on parse.
    position : int
        Genomic position (1-based, VCF convention).
    ref : str
        Reference allele.
    alt : str
        Alternate allele (first ALT for multi-allelic sites).
    sample_id : str
        Optional sample identifier from VCF header.
    genotype : str
        Genotype string (e.g. "0/1" for heterozygous).
    quality : float
        VCF QUAL score.
    consequence : Consequence
        Most severe VEP consequence for this variant.
    gene : str or None
        HGNC gene symbol from VEP.
    transcript : str or None
        Ensembl transcript ID from VEP.
    hgvsc : str or None
        HGVSc notation (coding sequence change).
    hgvsp : str or None
        HGVSp notation (protein change).
    sift : str or None
        SIFT prediction string (e.g. "deleterious").
    polyphen : str or None
        PolyPhen-2 prediction string (e.g. "probably_damaging").
    cadd_phred : float or None
        CADD PHRED score (higher = more deleterious).
    allele_frequency : float or None
        gnomAD v4.1 population allele frequency.
    gnomad_ac : int or None
        gnomAD allele count.
    gnomad_an : int or None
        gnomAD total alleles.
    clinvar_significance : str or None
        ClinVar clinical significance string.
    clinvar_id : str or None
        ClinVar variation ID.
    clinvar_conditions : list of str
        Associated disease conditions from ClinVar.
    annotation : dict
        Catch-all dict for any additional annotation fields.
    priority_score : float
        Composite priority score computed by VariantPrioritizer (0–100).
    """

    # ── Core VCF fields ───────────────────────────────────────────────────────
    chromosome: str
    position: int
    ref: str
    alt: str
    sample_id: str = "."
    genotype: str = "."
    quality: float = 0.0

    # ── VEP annotation ────────────────────────────────────────────────────────
    consequence: Consequence = Consequence.UNKNOWN
    gene: Optional[str] = None
    transcript: Optional[str] = None
    hgvsc: Optional[str] = None
    hgvsp: Optional[str] = None
    sift: Optional[str] = None
    polyphen: Optional[str] = None
    cadd_phred: Optional[float] = None

    # ── Population frequency (gnomAD v4.1) ───────────────────────────────────
    allele_frequency: Optional[float] = None
    gnomad_ac: Optional[int] = None
    gnomad_an: Optional[int] = None

    # ── ClinVar ───────────────────────────────────────────────────────────────
    clinvar_significance: Optional[str] = None
    clinvar_id: Optional[str] = None
    clinvar_conditions: List[str] = field(default_factory=list)

    # ── Catch-all and scoring ─────────────────────────────────────────────────
    annotation: Dict[str, Any] = field(default_factory=dict)
    priority_score: float = 0.0

    # ── Validation ────────────────────────────────────────────────────────────

    def __post_init__(self) -> None:
        if not self.chromosome:
            raise ValueError("Chromosome must be specified")
        if self.position <= 0:
            raise ValueError("Position must be positive")
        if not self.ref or not self.alt:
            raise ValueError("Reference and alternate alleles must be specified")

    # ── Derived properties ────────────────────────────────────────────────────

    @property
    def variant_id(self) -> str:
        """Canonical ID: chr-pos-ref-alt (matches VEP region notation output)."""
        return f"{self.chromosome}-{self.position}-{self.ref}-{self.alt}"

    def infer_type(self) -> VariantType:
        """Infer variant type from ref/alt lengths."""
        if len(self.ref) == len(self.alt) == 1:
            return VariantType.SNP
        elif len(self.ref) > len(self.alt):
            return VariantType.DELETION
        elif len(self.alt) > len(self.ref):
            return VariantType.INSERTION
        else:
            return VariantType.INDEL

    def is_rare(self, threshold: float = 0.01) -> bool:
        """
        Return True if the variant is rare relative to *threshold*.

        A variant with no known gnomAD frequency is treated as potentially
        rare (returns True), which is the conservative clinical assumption.
        """
        if self.allele_frequency is None:
            return True   # unknown frequency → assume rare
        return self.allele_frequency < threshold

    def is_potentially_pathogenic(self) -> bool:
        """Return True if the consequence is in the high-impact tier."""
        high_impact = {
            Consequence.TRANSCRIPT_ABLATION,
            Consequence.FRAMESHIFT,
            Consequence.STOP_GAINED,
            Consequence.STOP_LOST,
            Consequence.START_LOST,
            Consequence.SPLICE_ACCEPTOR,
            Consequence.SPLICE_DONOR,
            Consequence.DISRUPTIVE_INFRAME,
        }
        return self.consequence in high_impact

    def to_dict(self) -> Dict[str, Any]:
        """Serialise to a plain dict (enum values converted to strings)."""
        d = asdict(self)
        d["consequence"] = self.consequence.value if self.consequence else None
        return d

    def __str__(self) -> str:
        return f"{self.chromosome}:{self.position} {self.ref}>{self.alt}"

    def __repr__(self) -> str:
        return (
            f"Variant(chr={self.chromosome}, pos={self.position}, "
            f"{self.ref}>{self.alt}, gene={self.gene}, "
            f"consequence={self.consequence.value if self.consequence else 'unknown'})"
        )


@dataclass
class VariantCollection:
    """
    Container for managing multiple Variant objects.

    Supports filtering, sorting, and counting operations.
    All filter/sort methods return new VariantCollection instances
    (functional style — originals are not mutated).
    """

    variants: List[Variant] = field(default_factory=list)
    name: str = "variants"

    # ── Mutation ──────────────────────────────────────────────────────────────

    def add(self, variant: Variant) -> None:
        self.variants.append(variant)

    def add_many(self, variants: List[Variant]) -> None:
        self.variants.extend(variants)

    # ── Filters ───────────────────────────────────────────────────────────────

    def filter_by_consequence(self, consequence: Consequence) -> "VariantCollection":
        filtered = VariantCollection(name=f"{self.name}_filtered_{consequence.value}")
        filtered.add_many([v for v in self.variants if v.consequence == consequence])
        return filtered

    def filter_by_rarity(self, threshold: float = 0.01) -> "VariantCollection":
        filtered = VariantCollection(name=f"{self.name}_rare")
        filtered.add_many([v for v in self.variants if v.is_rare(threshold)])
        return filtered

    def filter_by_gene(self, gene: str) -> "VariantCollection":
        filtered = VariantCollection(name=f"{self.name}_{gene}")
        filtered.add_many([v for v in self.variants if v.gene == gene])
        return filtered

    # ── Sorting ───────────────────────────────────────────────────────────────

    def sort_by_priority(self) -> "VariantCollection":
        sorted_variants = sorted(self.variants, key=lambda v: v.priority_score, reverse=True)
        result = VariantCollection(name=f"{self.name}_sorted")
        result.add_many(sorted_variants)
        return result

    def top_n(self, n: int) -> "VariantCollection":
        top = sorted(self.variants, key=lambda v: v.priority_score, reverse=True)[:n]
        result = VariantCollection(name=f"{self.name}_top_{n}")
        result.add_many(top)
        return result

    # ── Statistics ────────────────────────────────────────────────────────────

    def count_by_gene(self) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for v in self.variants:
            gene = v.gene or "intergenic"
            counts[gene] = counts.get(gene, 0) + 1
        return counts

    def to_list(self) -> List[Dict[str, Any]]:
        return [v.to_dict() for v in self.variants]

    # ── Dunder ────────────────────────────────────────────────────────────────

    def __len__(self) -> int:
        return len(self.variants)

    def __iter__(self):
        return iter(self.variants)

    def __repr__(self) -> str:
        return f"VariantCollection(name={self.name!r}, variants={len(self.variants)})"

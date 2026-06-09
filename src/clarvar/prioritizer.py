"""
Variant prioritisation module.

Your task: score and rank annotated variants by predicted clinical importance.

Scoring model
-------------
Each variant receives a priority_score between 0 and 100 composed of three
weighted components:

    priority_score = (consequence_weight × consequence_score)
                   + (clinvar_weight     × clinvar_score)
                   + (frequency_weight   × frequency_score)

Default weights: consequence=0.4, clinvar=0.4, frequency=0.2
(Weights are normalised to sum to 1.0 regardless of input values.)

Component scoring tables
------------------------

Consequence score (0–100):
    frameshift_variant    → 100
    stop_gained           → 100
    splice_acceptor       → 95
    splice_donor          → 95
    stop_lost             → 95
    start_lost            → 90
    inframe_deletion      → 70
    inframe_insertion     → 70
    missense_variant      → 40
    synonymous_variant    → 10
    upstream/downstream   → 5
    intergenic            → 0
    unknown               → 20  (conservative default)

ClinVar score (0–100):
    pathogenic            → 100
    likely_pathogenic     → 80
    uncertain_significance → 30
    likely_benign         → 10
    benign                → 5
    no entry              → 20  (neutral default)

Frequency score (0–100):
    AF < 0.0001           → 100  (very rare)
    AF 0.0001–0.01        → 80   (rare)
    AF 0.01–0.05          → 60   (low frequency)
    AF 0.05–0.1           → 40
    AF > 0.1              → 20   (common)
    no AF data            → 50   (neutral — unknown is not benign)

Ranking strategies
------------------
Three pre-built strategies alter the weights for different use cases:
- ConsequenceOnlyStrategy  — consequence=1.0, others=0.0
- RarityFirstStrategy      — frequency=0.6, others=0.2 each
- ClinicalEvidenceFirstStrategy — clinvar=0.6, others=0.2 each
"""

import logging
from typing import Dict, List, Optional

from clarvar.variant import Consequence, Variant, VariantCollection

logger = logging.getLogger(__name__)


class VariantPrioritizer:
    """
    Scores and ranks variants by predicted clinical importance.

    Parameters
    ----------
    consequence_weight : float
        Weight for consequence severity (default 0.4).
    clinvar_weight : float
        Weight for ClinVar significance (default 0.4).
    frequency_weight : float
        Weight for population rarity (default 0.2).

    Notes
    -----
    Weights are normalised to sum to 1.0 on initialisation, so passing
    (0.5, 0.3, 0.2) is equivalent to (5, 3, 2).

    Example
    -------
    >>> prioritizer = VariantPrioritizer()
    >>> prioritized = prioritizer.prioritize_collection(annotated_variants)
    >>> report = prioritizer.generate_report(prioritized)
    """

    # ── Scoring tables — fill these in ───────────────────────────────────────

    CONSEQUENCE_SCORES: Dict[Consequence, float] = {
        # TODO: populate from the scoring model in the module docstring
    }

    CLINVAR_SCORES: Dict[str, float] = {
        # TODO: populate from the scoring model in the module docstring
        # Keys should be lowercase ClinVar significance strings
    }

    def __init__(
        self,
        consequence_weight: float = 0.4,
        clinvar_weight: float = 0.4,
        frequency_weight: float = 0.2,
    ):
        total = consequence_weight + clinvar_weight + frequency_weight
        self.consequence_weight = consequence_weight / total
        self.clinvar_weight = clinvar_weight / total
        self.frequency_weight = frequency_weight / total

        logger.info(
            f"Initialized prioritizer with weights: "
            f"consequence={self.consequence_weight:.2f}, "
            f"clinvar={self.clinvar_weight:.2f}, "
            f"frequency={self.frequency_weight:.2f}"
        )

    # ── Component scoring — implement these ───────────────────────────────────

    def _score_consequence(self, variant: Variant) -> float:
        """
        Return a consequence severity score (0–100) for the variant.

        Look up variant.consequence in CONSEQUENCE_SCORES.
        Return 20 as a conservative default for UNKNOWN or missing keys.

        Parameters
        ----------
        variant : Variant

        Returns
        -------
        float
        """
        raise NotImplementedError("TODO: implement _score_consequence")

    def _score_clinvar(self, variant: Variant) -> float:
        """
        Return a ClinVar significance score (0–100) for the variant.

        Look up variant.clinvar_significance (lowercased) in CLINVAR_SCORES.
        Return 20 as a neutral default if no ClinVar entry is present.

        Parameters
        ----------
        variant : Variant

        Returns
        -------
        float
        """
        raise NotImplementedError("TODO: implement _score_clinvar")

    def _score_frequency(self, variant: Variant) -> float:
        """
        Return a population rarity score (0–100) for the variant.

        Use the frequency tiers defined in the module docstring.
        Return 50 (neutral) when variant.allele_frequency is None —
        unknown frequency is not the same as common.

        Parameters
        ----------
        variant : Variant

        Returns
        -------
        float
        """
        raise NotImplementedError("TODO: implement _score_frequency")

    # ── Main prioritisation methods ───────────────────────────────────────────

    def prioritize_variant(self, variant: Variant) -> Variant:
        """
        Compute and set the priority_score on a single Variant.

        priority_score = (consequence_weight × _score_consequence())
                       + (clinvar_weight     × _score_clinvar())
                       + (frequency_weight   × _score_frequency())

        Store the result in variant.priority_score and return the variant.

        Parameters
        ----------
        variant : Variant

        Returns
        -------
        Variant
            Same object with priority_score set.
        """
        raise NotImplementedError("TODO: implement prioritize_variant")

    def prioritize_collection(
        self, collection: VariantCollection
    ) -> VariantCollection:
        """
        Score all variants and return them sorted by priority (highest first).

        Parameters
        ----------
        collection : VariantCollection

        Returns
        -------
        VariantCollection
            New collection with all variants scored and sorted descending
            by priority_score.
        """
        raise NotImplementedError("TODO: implement prioritize_collection")

    # ── Reporting ─────────────────────────────────────────────────────────────

    def generate_report(self, collection: VariantCollection) -> Dict:
        """
        Generate a summary dict of prioritisation results.

        Returns a dict with the following keys:
            total_variants       : int
            mean_priority        : float
            median_priority      : float
            max_priority         : float
            min_priority         : float
            top_variants         : list of dicts (up to 20), each with keys:
                                   variant, gene, consequence, clinvar,
                                   frequency, score
            consequence_distribution : dict mapping consequence string → count
            clinvar_distribution     : dict mapping significance string → count

        Return {"total_variants": 0, "mean_priority": 0, "top_variants": []}
        if the collection is empty.

        Parameters
        ----------
        collection : VariantCollection

        Returns
        -------
        dict
        """
        raise NotImplementedError("TODO: implement generate_report")

    def _get_consequence_distribution(
        self, collection: VariantCollection
    ) -> Dict[str, int]:
        """Return a count of variants by consequence term string."""
        raise NotImplementedError("TODO: implement _get_consequence_distribution")

    def _get_clinvar_distribution(
        self, collection: VariantCollection
    ) -> Dict[str, int]:
        """Return a count of variants by ClinVar significance string."""
        raise NotImplementedError("TODO: implement _get_clinvar_distribution")


# ── Ranking strategies ────────────────────────────────────────────────────────

class RankingStrategy:
    """Base class for alternative ranking strategies."""

    def rank(self, variants: List[Variant]) -> List[Variant]:
        """
        Rank a list of variants and return them sorted highest-first.

        Parameters
        ----------
        variants : list of Variant

        Returns
        -------
        list of Variant
        """
        raise NotImplementedError("Subclasses must implement rank()")


class ConsequenceOnlyStrategy(RankingStrategy):
    """
    Ranks variants by functional consequence severity alone.

    Weights: consequence=1.0, clinvar=0.0, frequency=0.0
    Use this when you want a purely functional impact ranking with no
    influence from population frequency or clinical databases.
    """

    def rank(self, variants: List[Variant]) -> List[Variant]:
        raise NotImplementedError("TODO: implement ConsequenceOnlyStrategy.rank")


class RarityFirstStrategy(RankingStrategy):
    """
    Prioritises rare variants above all others.

    Weights: consequence=0.2, clinvar=0.2, frequency=0.6
    Use this for rare disease workflows where population rarity is the
    primary filter.
    """

    def rank(self, variants: List[Variant]) -> List[Variant]:
        raise NotImplementedError("TODO: implement RarityFirstStrategy.rank")


class ClinicalEvidenceFirstStrategy(RankingStrategy):
    """
    Prioritises variants with strong ClinVar clinical evidence.

    Weights: consequence=0.2, clinvar=0.6, frequency=0.2
    Use this when you want to surface previously classified pathogenic
    variants regardless of population frequency.
    """

    def rank(self, variants: List[Variant]) -> List[Variant]:
        raise NotImplementedError(
            "TODO: implement ClinicalEvidenceFirstStrategy.rank"
        )

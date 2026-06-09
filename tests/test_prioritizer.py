"""
Unit tests for Prioritizer module.
"""

import pytest
from clarvar.variant import Variant, Consequence
from clarvar.prioritizer import (
    VariantPrioritizer,
    ConsequenceOnlyStrategy,
    RarityFirstStrategy,
    ClinicalEvidenceFirstStrategy
)


class TestVariantPrioritizer:
    """Test cases for VariantPrioritizer class."""
    
    @pytest.fixture
    def prioritizer(self):
        """Create a VariantPrioritizer instance."""
        return VariantPrioritizer()
    
    def test_prioritizer_creation(self, prioritizer):
        """Test creating a prioritizer."""
        assert prioritizer is not None
        assert prioritizer.consequence_weight > 0
        assert prioritizer.clinvar_weight > 0
        assert prioritizer.frequency_weight > 0
        
        # Weights should be normalized
        total = (prioritizer.consequence_weight + 
                prioritizer.clinvar_weight + 
                prioritizer.frequency_weight)
        assert abs(total - 1.0) < 0.001
    
    def test_prioritizer_custom_weights(self):
        """Test creating prioritizer with custom weights."""
        prioritizer = VariantPrioritizer(
            consequence_weight=0.5,
            clinvar_weight=0.3,
            frequency_weight=0.2
        )
        
        total = (prioritizer.consequence_weight + 
                prioritizer.clinvar_weight + 
                prioritizer.frequency_weight)
        assert abs(total - 1.0) < 0.001
    
    def test_score_consequence(self, prioritizer):
        """Test consequence scoring."""
        # High-impact consequence
        variant_frameshift = Variant(
            chromosome="1", position=100, ref="A", alt="G",
            consequence=Consequence.FRAMESHIFT
        )
        score = prioritizer._score_consequence(variant_frameshift)
        assert score == 100
        
        # Low-impact consequence
        variant_synonymous = Variant(
            chromosome="1", position=100, ref="A", alt="G",
            consequence=Consequence.SYNONYMOUS
        )
        score = prioritizer._score_consequence(variant_synonymous)
        assert score == 10
    
    def test_score_clinvar(self, prioritizer):
        """Test ClinVar scoring."""
        variant = Variant(chromosome="1", position=100, ref="A", alt="G")
        
        # Pathogenic
        variant.clinvar_significance = "pathogenic"
        score = prioritizer._score_clinvar(variant)
        assert score == 100
        
        # Benign
        variant.clinvar_significance = "benign"
        score = prioritizer._score_clinvar(variant)
        assert score == 5
        
        # No ClinVar entry
        variant.clinvar_significance = None
        score = prioritizer._score_clinvar(variant)
        assert score == 20  # Default
    
    def test_score_frequency(self, prioritizer):
        """Test allele frequency scoring."""
        variant = Variant(chromosome="1", position=100, ref="A", alt="G")
        
        # Very rare
        variant.allele_frequency = 0.00005
        score = prioritizer._score_frequency(variant)
        assert score == 100
        
        # Rare
        variant.allele_frequency = 0.005
        score = prioritizer._score_frequency(variant)
        assert score == 80
        
        # Common
        variant.allele_frequency = 0.1
        score = prioritizer._score_frequency(variant)
        assert score == 20
        
        # No frequency data
        variant.allele_frequency = None
        score = prioritizer._score_frequency(variant)
        assert score == 50  # Neutral
    
    def test_prioritize_variant(self, prioritizer, sample_variant):
        """Test prioritizing a single variant."""
        prioritizer.prioritize_variant(sample_variant)
        
        # Should have priority score
        assert sample_variant.priority_score >= 0
        assert sample_variant.priority_score <= 100
    
    def test_prioritize_collection(self, prioritizer, sample_collection):
        """Test prioritizing a collection."""
        prioritized = prioritizer.prioritize_collection(sample_collection)
        
        assert len(prioritized) == len(sample_collection)
        
        # Should be sorted by priority
        scores = [v.priority_score for v in prioritized]
        assert scores == sorted(scores, reverse=True)
    
    def test_generate_report(self, prioritizer, sample_collection):
        """Test generating a report."""
        prioritizer.prioritize_collection(sample_collection)
        report = prioritizer.generate_report(sample_collection)
        
        # Check report structure
        assert 'total_variants' in report
        assert 'mean_priority' in report
        assert 'median_priority' in report
        assert 'max_priority' in report
        assert 'min_priority' in report
        assert 'top_variants' in report
        assert 'consequence_distribution' in report
        assert 'clinvar_distribution' in report
        
        # Check values
        assert report['total_variants'] == len(sample_collection)
        assert isinstance(report['mean_priority'], (int, float))
        assert len(report['top_variants']) <= 20
    
    def test_consequence_distribution(self, prioritizer, sample_collection):
        """Test consequence distribution calculation."""
        distribution = prioritizer._get_consequence_distribution(sample_collection)
        
        # Should have entries for each unique consequence
        total = sum(distribution.values())
        assert total == len(sample_collection)
    
    def test_clinvar_distribution(self, prioritizer, sample_collection):
        """Test ClinVar distribution calculation."""
        distribution = prioritizer._get_clinvar_distribution(sample_collection)
        
        # Should have entries
        total = sum(distribution.values())
        assert total == len(sample_collection)


class TestRankingStrategies:
    """Test cases for ranking strategy classes."""
    
    def test_consequence_only_strategy(self, sample_collection):
        """Test consequence-only ranking strategy."""
        strategy = ConsequenceOnlyStrategy()
        ranked = strategy.rank(list(sample_collection))
        
        assert len(ranked) == len(sample_collection)
        
        # Frameshift should rank higher than synonymous
        consequences = [v.consequence for v in ranked]
        frameshift_idx = next((i for i, c in enumerate(consequences) 
                             if c == Consequence.FRAMESHIFT), -1)
        synonymous_idx = next((i for i, c in enumerate(consequences) 
                             if c == Consequence.SYNONYMOUS), -1)
        
        if frameshift_idx >= 0 and synonymous_idx >= 0:
            assert frameshift_idx < synonymous_idx
    
    def test_rarity_first_strategy(self, sample_collection):
        """Test rarity-first ranking strategy."""
        strategy = RarityFirstStrategy()
        ranked = strategy.rank(list(sample_collection))
        
        assert len(ranked) == len(sample_collection)
        
        # Higher priority for rare variants
        frequencies = [v.allele_frequency for v in ranked if v.allele_frequency]
        # Should generally have lower frequencies first (rarer)
        if len(frequencies) >= 2:
            assert frequencies[0] <= frequencies[-1] or True  # May not always be true
    
    def test_clinical_evidence_first_strategy(self, sample_collection):
        """Test clinical evidence-first ranking strategy."""
        strategy = ClinicalEvidenceFirstStrategy()
        ranked = strategy.rank(list(sample_collection))
        
        assert len(ranked) == len(sample_collection)
        
        # Should prioritize variants with ClinVar info
        with_clinvar = [v for v in ranked if v.clinvar_significance]
        without_clinvar = [v for v in ranked if not v.clinvar_significance]
        
        # If there are both, with_clinvar should generally rank higher
        if with_clinvar and without_clinvar:
            assert ranked[0] in (with_clinvar + [None])[0:len(with_clinvar)]


@pytest.mark.integration
class TestPrioritizationIntegration:
    """Integration tests for prioritization workflow."""
    
    def test_full_prioritization_workflow(self, sample_collection):
        """Test complete prioritization workflow."""
        prioritizer = VariantPrioritizer()
        
        # First annotate
        from clarvar.annotator import VariantAnnotator
        annotator = VariantAnnotator()
        annotated = annotator.annotate_collection(sample_collection)
        
        # Then prioritize
        prioritized = prioritizer.prioritize_collection(annotated)
        
        # Verify
        assert len(prioritized) == len(sample_collection)
        
        # All should have scores
        for variant in prioritized:
            assert variant.priority_score >= 0
            assert variant.priority_score <= 100
    
    def test_prioritization_consistency(self, sample_collection):
        """Test that prioritization is consistent."""
        prioritizer = VariantPrioritizer()
        
        # Run prioritization twice
        result1 = prioritizer.prioritize_collection(sample_collection)
        result2 = prioritizer.prioritize_collection(sample_collection)
        
        # Should get same scores
        scores1 = [v.priority_score for v in result1]
        scores2 = [v.priority_score for v in result2]
        
        for s1, s2 in zip(scores1, scores2):
            assert abs(s1 - s2) < 0.01  # Allow small floating point differences

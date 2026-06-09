"""
Unit tests for Variant class.
"""

import pytest
from clarvar.variant import Variant, VariantCollection, VariantType, Consequence


class TestVariant:
    """Test cases for Variant class."""
    
    def test_variant_creation(self):
        """Test basic variant creation."""
        variant = Variant(
            chromosome="1",
            position=12345,
            ref="A",
            alt="G",
            gene="TP53"
        )
        
        assert variant.chromosome == "1"
        assert variant.position == 12345
        assert variant.ref == "A"
        assert variant.alt == "G"
        assert variant.gene == "TP53"
    
    def test_variant_validation(self):
        """Test variant validation."""
        # Invalid position
        with pytest.raises(ValueError):
            Variant(chromosome="1", position=0, ref="A", alt="G")
        
        # Missing chromosome
        with pytest.raises(ValueError):
            Variant(chromosome="", position=100, ref="A", alt="G")
        
        # Missing alleles
        with pytest.raises(ValueError):
            Variant(chromosome="1", position=100, ref="", alt="G")
    
    def test_infer_type(self):
        """Test variant type inference."""
        # SNP
        snp = Variant(chromosome="1", position=100, ref="A", alt="G")
        assert snp.infer_type() == VariantType.SNP
        
        # Deletion
        deletion = Variant(chromosome="1", position=100, ref="ATG", alt="A")
        assert deletion.infer_type() == VariantType.DELETION
        
        # Insertion
        insertion = Variant(chromosome="1", position=100, ref="A", alt="ATG")
        assert insertion.infer_type() == VariantType.INSERTION
    
    def test_is_rare(self):
        """Test rarity assessment."""
        variant = Variant(chromosome="1", position=100, ref="A", alt="G")
        
        # No frequency data — clinically conservative: treat as potentially rare
        assert variant.is_rare()
        
        # Rare variant
        variant.allele_frequency = 0.0001
        assert variant.is_rare(threshold=0.01)
        
        # Common variant
        variant.allele_frequency = 0.1
        assert not variant.is_rare(threshold=0.01)
    
    def test_is_potentially_pathogenic(self):
        """Test pathogenicity assessment."""
        variant = Variant(chromosome="1", position=100, ref="A", alt="G")
        
        # Frameshift
        variant.consequence = Consequence.FRAMESHIFT
        assert variant.is_potentially_pathogenic()
        
        # Missense (not pathogenic by consequence alone)
        variant.consequence = Consequence.MISSENSE
        assert not variant.is_potentially_pathogenic()
        
        # Splice donor
        variant.consequence = Consequence.SPLICE_DONOR
        assert variant.is_potentially_pathogenic()
    
    def test_variant_to_dict(self):
        """Test variant serialization."""
        variant = Variant(
            chromosome="17",
            position=41246881,
            ref="G",
            alt="A",
            gene="BRCA1",
            consequence=Consequence.STOP_GAINED,
            priority_score=95.5
        )
        
        d = variant.to_dict()
        assert d["chromosome"] == "17"
        assert d["position"] == 41246881
        assert d["gene"] == "BRCA1"
        assert d["consequence"] == "stop_gained"
        assert d["priority_score"] == 95.5
    
    def test_variant_string_representation(self):
        """Test variant string representations."""
        variant = Variant(chromosome="1", position=100, ref="A", alt="G")
        
        assert str(variant) == "1:100 A>G"
        assert "chr=" in repr(variant)


class TestVariantCollection:
    """Test cases for VariantCollection class."""
    
    @pytest.fixture
    def sample_variants(self):
        """Create sample variants for testing."""
        variants = [
            Variant(chromosome="1", position=100, ref="A", alt="G", gene="GENE1",
                   consequence=Consequence.MISSENSE, allele_frequency=0.05),
            Variant(chromosome="1", position=200, ref="T", alt="C", gene="GENE2",
                   consequence=Consequence.FRAMESHIFT, allele_frequency=0.001),
            Variant(chromosome="2", position=300, ref="G", alt="A", gene="GENE1",
                   consequence=Consequence.SYNONYMOUS, allele_frequency=0.2),
        ]
        return variants
    
    def test_collection_creation(self):
        """Test collection creation."""
        collection = VariantCollection(name="test_variants")
        assert len(collection) == 0
        assert collection.name == "test_variants"
    
    def test_add_variant(self, sample_variants):
        """Test adding variants."""
        collection = VariantCollection()
        collection.add(sample_variants[0])
        
        assert len(collection) == 1
        assert collection.variants[0].gene == "GENE1"
    
    def test_add_many(self, sample_variants):
        """Test adding multiple variants."""
        collection = VariantCollection()
        collection.add_many(sample_variants)
        
        assert len(collection) == 3
    
    def test_filter_by_consequence(self, sample_variants):
        """Test filtering by consequence."""
        collection = VariantCollection()
        collection.add_many(sample_variants)
        
        frameshift_only = collection.filter_by_consequence(Consequence.FRAMESHIFT)
        assert len(frameshift_only) == 1
        assert frameshift_only.variants[0].position == 200
    
    def test_filter_by_rarity(self, sample_variants):
        """Test filtering by rarity."""
        collection = VariantCollection()
        collection.add_many(sample_variants)
        
        rare = collection.filter_by_rarity(threshold=0.01)
        assert len(rare) == 1
        assert rare.variants[0].position == 200
    
    def test_filter_by_gene(self, sample_variants):
        """Test filtering by gene."""
        collection = VariantCollection()
        collection.add_many(sample_variants)
        
        gene1_variants = collection.filter_by_gene("GENE1")
        assert len(gene1_variants) == 2
    
    def test_sort_by_priority(self, sample_variants):
        """Test sorting by priority."""
        # Set priority scores
        sample_variants[0].priority_score = 10
        sample_variants[1].priority_score = 90
        sample_variants[2].priority_score = 50
        
        collection = VariantCollection()
        collection.add_many(sample_variants)
        
        sorted_collection = collection.sort_by_priority()
        assert sorted_collection.variants[0].priority_score == 90
        assert sorted_collection.variants[2].priority_score == 10
    
    def test_top_n(self, sample_variants):
        """Test getting top N variants."""
        # Set priority scores
        for i, variant in enumerate(sample_variants):
            variant.priority_score = (i + 1) * 10
        
        collection = VariantCollection()
        collection.add_many(sample_variants)
        
        top2 = collection.top_n(2)
        assert len(top2) == 2
        assert top2.variants[0].priority_score == 30
    
    def test_count_by_gene(self, sample_variants):
        """Test counting variants by gene."""
        collection = VariantCollection()
        collection.add_many(sample_variants)
        
        counts = collection.count_by_gene()
        assert counts["GENE1"] == 2
        assert counts["GENE2"] == 1
    
    def test_collection_iteration(self, sample_variants):
        """Test iterating over collection."""
        collection = VariantCollection()
        collection.add_many(sample_variants)
        
        count = 0
        for variant in collection:
            assert isinstance(variant, Variant)
            count += 1
        
        assert count == 3


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

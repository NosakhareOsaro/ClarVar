"""
Pytest configuration and fixtures for ClarVar tests.
"""

import pytest
from pathlib import Path
from clarvar.variant import Variant, VariantCollection, Consequence


@pytest.fixture
def sample_variant():
    """Create a sample variant for testing."""
    return Variant(
        chromosome="1",
        position=12345,
        ref="A",
        alt="G",
        gene="TP53",
        consequence=Consequence.MISSENSE,
        allele_frequency=0.0001,
        clinvar_significance="uncertain_significance"
    )


@pytest.fixture
def sample_variants():
    """Create multiple sample variants for testing."""
    return [
        Variant(
            chromosome="1",
            position=100,
            ref="A",
            alt="G",
            gene="GENE1",
            consequence=Consequence.MISSENSE,
            allele_frequency=0.05
        ),
        Variant(
            chromosome="1",
            position=200,
            ref="T",
            alt="C",
            gene="GENE2",
            consequence=Consequence.FRAMESHIFT,
            allele_frequency=0.001
        ),
        Variant(
            chromosome="2",
            position=300,
            ref="G",
            alt="A",
            gene="GENE1",
            consequence=Consequence.SYNONYMOUS,
            allele_frequency=0.2
        ),
        Variant(
            chromosome="17",
            position=41246881,
            ref="G",
            alt="A",
            gene="BRCA1",
            consequence=Consequence.STOP_GAINED,
            allele_frequency=0.00001,
            clinvar_significance="pathogenic"
        ),
    ]


@pytest.fixture
def sample_collection(sample_variants):
    """Create a sample variant collection."""
    collection = VariantCollection(name="test_collection")
    collection.add_many(sample_variants)
    return collection


@pytest.fixture
def example_vcf_file():
    """Path to example VCF file."""
    vcf_path = Path(__file__).parent.parent / "examples" / "sample_variants.vcf"
    if vcf_path.exists():
        return vcf_path
    return None


@pytest.fixture
def temp_vcf_file(tmp_path):
    """Create a temporary VCF file for testing."""
    vcf_file = tmp_path / "test_variants.vcf"
    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=CSQ,Description="Consequence">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	100	PASS	.
1	200	.	T	C	95	PASS	.
2	300	.	G	A	88	PASS	.
"""
    vcf_file.write_text(vcf_content)
    return vcf_file


# Configuration
def pytest_configure(config):
    """Configure pytest."""
    config.addinivalue_line(
        "markers", "integration: mark test as an integration test"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow running"
    )


# Logging
@pytest.fixture(scope="session", autouse=True)
def setup_logging():
    """Set up logging for tests."""
    import logging
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Suppress verbose library logging during tests
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    logging.getLogger("requests").setLevel(logging.WARNING)

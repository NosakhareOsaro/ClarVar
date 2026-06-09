"""
Tests for the annotator module.

The VariantAnnotator now calls the real Ensembl VEP REST API, so
live-API tests are skipped in CI unless CLARVAR_LIVE_TESTS=1 is set.
The LocalAnnotator and all unit-level annotator logic are fully tested offline.
"""

import os
import pytest
from unittest.mock import patch, MagicMock

from clarvar.variant import Variant, Consequence
from clarvar.annotator import (
    VariantAnnotator, LocalAnnotator,
    _to_vep_region, _extract_gnomad_af, _extract_clinvar, _apply_vep_result,
)

LIVE = os.environ.get("CLARVAR_LIVE_TESTS") == "1"


# ── Helper ─────────────────────────────────────────────────────────────────────

def make_variant(**kwargs) -> Variant:
    defaults = dict(chromosome="1", position=100, ref="A", alt="G")
    defaults.update(kwargs)
    return Variant(**defaults)


# ── VEP region notation ────────────────────────────────────────────────────────

class TestVepRegionNotation:
    def test_snv(self):
        v = make_variant(chromosome="1", position=100, ref="A", alt="T")
        assert _to_vep_region(v) == "1 100 100 A/T 1"

    def test_deletion_end_pos(self):
        # ref=CTT alt=C → end = pos + len(ref) - 1 = 1000 + 3 - 1 = 1002
        v = make_variant(chromosome="7", position=1000, ref="CTT", alt="C")
        assert _to_vep_region(v) == "7 1000 1002 CTT/C 1"

    def test_insertion_same_end(self):
        v = make_variant(chromosome="2", position=500, ref="A", alt="ATG")
        assert _to_vep_region(v) == "2 500 500 A/ATG 1"


# ── gnomAD AF extraction ───────────────────────────────────────────────────────

class TestGnomadExtraction:
    def _hit(self, key, value):
        return {"colocated_variants": [{"frequencies": {"A": {key: value}}}]}

    def test_gnomad_key(self):
        assert _extract_gnomad_af(self._hit("gnomad", 0.001)) == pytest.approx(0.001)

    def test_gnomadg_key(self):
        assert _extract_gnomad_af(self._hit("gnomadg", 0.005)) == pytest.approx(0.005)

    def test_gnomade_key(self):
        assert _extract_gnomad_af(self._hit("gnomade", 0.0001)) == pytest.approx(0.0001)

    def test_no_colocated(self):
        assert _extract_gnomad_af({}) is None

    def test_no_gnomad_key(self):
        hit = {"colocated_variants": [{"frequencies": {"A": {"other": 0.1}}}]}
        assert _extract_gnomad_af(hit) is None


# ── ClinVar extraction ─────────────────────────────────────────────────────────

class TestClinvarExtraction:
    def test_pathogenic(self):
        hit = {
            "colocated_variants": [{
                "clin_sig": ["pathogenic"],
                "var_synonyms": {"ClinVar": ["12345"]},
            }]
        }
        sig, cid = _extract_clinvar(hit)
        assert sig == "pathogenic"
        assert cid == "12345"

    def test_list_significance(self):
        hit = {
            "colocated_variants": [{
                "clin_sig": ["likely_pathogenic", "pathogenic"],
            }]
        }
        sig, cid = _extract_clinvar(hit)
        assert "likely_pathogenic" in sig
        assert cid is None

    def test_somatic_skipped(self):
        hit = {
            "colocated_variants": [
                {"somatic": True, "clin_sig": ["pathogenic"]},
                {"clin_sig": ["benign"]},
            ]
        }
        sig, _ = _extract_clinvar(hit)
        assert sig == "benign"

    def test_no_clinvar(self):
        sig, cid = _extract_clinvar({})
        assert sig is None
        assert cid is None


# ── apply_vep_result ───────────────────────────────────────────────────────────

class TestApplyVepResult:
    def _mock_hit(self):
        return {
            "input": "1 100 100 A/T 1",
            "transcript_consequences": [{
                "gene_symbol": "TP53",
                "transcript_id": "ENST00000269305",
                "consequence_terms": ["missense_variant"],
                "hgvsc": "ENST00000269305.9:c.100A>T",
                "hgvsp": "ENSP00000269305.4:p.Arg34Trp",
                "sift_prediction": "deleterious",
                "polyphen_prediction": "probably_damaging",
                "cadd_phred": 28.5,
            }],
            "colocated_variants": [{
                "clin_sig": ["pathogenic"],
                "frequencies": {"T": {"gnomad": 0.0001}},
            }],
        }

    def test_gene_and_transcript(self):
        v = make_variant()
        _apply_vep_result(v, self._mock_hit())
        assert v.gene == "TP53"
        assert v.transcript == "ENST00000269305"

    def test_consequence_parsed(self):
        v = make_variant()
        _apply_vep_result(v, self._mock_hit())
        assert v.consequence == Consequence.MISSENSE

    def test_hgvs_fields(self):
        v = make_variant()
        _apply_vep_result(v, self._mock_hit())
        assert "c.100A>T" in v.hgvsc
        assert "p.Arg34Trp" in v.hgvsp

    def test_sift_polyphen(self):
        v = make_variant()
        _apply_vep_result(v, self._mock_hit())
        assert v.sift == "deleterious"
        assert v.polyphen == "probably_damaging"

    def test_cadd(self):
        v = make_variant()
        _apply_vep_result(v, self._mock_hit())
        assert v.cadd_phred == pytest.approx(28.5)

    def test_gnomad_af(self):
        v = make_variant()
        _apply_vep_result(v, self._mock_hit())
        assert v.allele_frequency == pytest.approx(0.0001)

    def test_clinvar_sig(self):
        v = make_variant()
        _apply_vep_result(v, self._mock_hit())
        assert v.clinvar_significance == "pathogenic"


# ── VariantAnnotator (mocked API) ─────────────────────────────────────────────

class TestVariantAnnotatorMocked:
    """Test the annotator class with the VEP API mocked out."""

    @patch("clarvar.annotator._post_vep_batch")
    def test_annotate_variant_calls_vep(self, mock_post):
        mock_post.return_value = [{
            "input": "1 100 100 A/G 1",
            "transcript_consequences": [{
                "gene_symbol": "BRCA1",
                "consequence_terms": ["missense_variant"],
            }],
            "colocated_variants": [],
        }]
        v = make_variant()
        annotator = VariantAnnotator()
        result = annotator.annotate_variant(v)
        assert result.gene == "BRCA1"
        assert result.consequence == Consequence.MISSENSE
        mock_post.assert_called_once()

    @patch("clarvar.annotator._post_vep_batch")
    def test_annotate_collection(self, mock_post, sample_collection):
        mock_post.return_value = []  # empty → no annotation, but no crash
        annotator = VariantAnnotator()
        result = annotator.annotate_collection(sample_collection)
        assert len(result) == len(sample_collection)

    def test_clear_cache(self):
        annotator = VariantAnnotator()
        annotator._response_cache["key"] = "value"
        annotator.clear_cache()
        assert len(annotator._response_cache) == 0

    def test_grch37_endpoint(self):
        annotator = VariantAnnotator(assembly="GRCh37")
        assert "grch37" in annotator._endpoint

    @patch("clarvar.annotator._post_vep_batch")
    def test_annotation_preserves_core_fields(self, mock_post, sample_variant):
        mock_post.return_value = []
        original_chrom = sample_variant.chromosome
        original_pos = sample_variant.position
        original_ref = sample_variant.ref
        original_alt = sample_variant.alt

        annotator = VariantAnnotator()
        annotator.annotate_variant(sample_variant)

        assert sample_variant.chromosome == original_chrom
        assert sample_variant.position == original_pos
        assert sample_variant.ref == original_ref
        assert sample_variant.alt == original_alt


# ── LocalAnnotator ────────────────────────────────────────────────────────────

class TestLocalAnnotator:
    def test_creation(self):
        la = LocalAnnotator()
        assert la is not None

    def test_load_databases_no_crash(self, tmp_path):
        la = LocalAnnotator(db_path=str(tmp_path / "db"))
        la.load_databases()  # should not raise

    def test_annotate_with_empty_db(self, sample_variant):
        la = LocalAnnotator()
        result = la.annotate_variant(sample_variant)
        assert result is sample_variant  # same object returned

    def test_annotate_with_gnomad_data(self, sample_variant):
        la = LocalAnnotator()
        key = f"{sample_variant.chromosome}:{sample_variant.position}_{sample_variant.ref}_{sample_variant.alt}"
        la.gnomad_db[key] = {"af": 0.0005}
        la.annotate_variant(sample_variant)
        assert sample_variant.allele_frequency == pytest.approx(0.0005)

    def test_annotate_with_clinvar_data(self, sample_variant):
        la = LocalAnnotator()
        key = f"{sample_variant.chromosome}:{sample_variant.position}_{sample_variant.ref}_{sample_variant.alt}"
        la.clinvar_db[key] = {"significance": "pathogenic"}
        la.annotate_variant(sample_variant)
        assert sample_variant.clinvar_significance == "pathogenic"

    def test_annotate_collection(self, sample_collection):
        la = LocalAnnotator()
        result = la.annotate_collection(sample_collection)
        assert len(result) == len(sample_collection)


# ── Live API integration test (skipped unless env var set) ────────────────────

@pytest.mark.skipif(not LIVE, reason="Set CLARVAR_LIVE_TESTS=1 to run live API tests")
class TestAnnotationLive:
    def test_annotate_known_pathogenic_variant(self):
        """Annotate BRCA2 c.5946delT (NM_000059.3), a known pathogenic frameshift."""
        v = Variant(chromosome="13", position=32914437, ref="AT", alt="A")
        annotator = VariantAnnotator(assembly="GRCh38")
        annotator.annotate_variant(v)

        assert v.gene is not None, "Gene should be populated"
        assert v.consequence != Consequence.UNKNOWN, "Consequence should be populated"
        # BRCA2 frameshift should have high CADD
        if v.cadd_phred is not None:
            assert v.cadd_phred > 20

"""
Tests for the reporter module.

Run with: pytest tests/test_reporter.py -v

Your task: make all these tests pass by implementing write_tsv(),
write_html_report(), _fmt(), _clinvar_badge(), and _priority_bar()
in src/clarvar/reporter.py
"""

import csv
from pathlib import Path

import pytest

from clarvar.variant import Variant, Consequence, VariantCollection
from clarvar.reporter import (
    write_tsv,
    write_html_report,
    _fmt,
    _clinvar_badge,
    _priority_bar,
    TSV_COLUMNS,
)


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def sample_variant():
    v = Variant(chromosome="17", position=41246881, ref="G", alt="A", gene="BRCA1")
    v.consequence = Consequence.MISSENSE
    v.hgvsc = "ENST00000357654.9:c.1067G>A"
    v.hgvsp = "ENSP00000350283.3:p.Arg356His"
    v.allele_frequency = 0.0001
    v.cadd_phred = 28.5
    v.clinvar_significance = "pathogenic"
    v.clinvar_id = "12345"
    v.sift = "deleterious"
    v.polyphen = "probably_damaging"
    v.priority_score = 85.0
    return v


@pytest.fixture
def sample_collection(sample_variant):
    v2 = Variant(chromosome="1", position=100, ref="A", alt="T")
    v2.consequence = Consequence.SYNONYMOUS
    v2.allele_frequency = 0.15
    v2.priority_score = 12.0
    col = VariantCollection()
    col.add(sample_variant)
    col.add(v2)
    return col


# ── _fmt ─────────────────────────────────────────────────────────────────────

class TestFmt:
    def test_none_returns_dot(self):
        assert _fmt(None) == "."

    def test_float_formatted(self):
        assert _fmt(0.123456, precision=2) == "0.12"

    def test_string_passthrough(self):
        assert _fmt("missense_variant") == "missense_variant"

    def test_int_as_string(self):
        assert _fmt(42) == "42"

    def test_zero_float(self):
        assert _fmt(0.0, precision=4) == "0.0000"


# ── _clinvar_badge ────────────────────────────────────────────────────────────

class TestClinvarBadge:
    def test_none_returns_empty(self):
        assert _clinvar_badge(None) == ""

    def test_pathogenic_is_red(self):
        badge = _clinvar_badge("pathogenic")
        assert "#c0392b" in badge
        assert "pathogenic" in badge

    def test_likely_pathogenic_is_orange(self):
        badge = _clinvar_badge("likely_pathogenic")
        assert "#e67e22" in badge

    def test_uncertain_is_purple(self):
        badge = _clinvar_badge("uncertain_significance")
        assert "#8e44ad" in badge

    def test_benign_is_green(self):
        assert "#27ae60" in _clinvar_badge("benign")
        assert "#27ae60" in _clinvar_badge("likely_benign")

    def test_unknown_sig_is_grey(self):
        badge = _clinvar_badge("conflicting_interpretations")
        assert "#7f8c8d" in badge

    def test_html_escaped(self):
        badge = _clinvar_badge("<script>")
        assert "<script>" not in badge


# ── _priority_bar ─────────────────────────────────────────────────────────────

class TestPriorityBar:
    def test_returns_html_string(self):
        bar = _priority_bar(50.0)
        assert "<div" in bar

    def test_high_score_is_red(self):
        assert "#e74c3c" in _priority_bar(75.0)

    def test_medium_score_is_orange(self):
        assert "#e67e22" in _priority_bar(45.0)

    def test_low_score_is_blue(self):
        assert "#3498db" in _priority_bar(20.0)

    def test_score_displayed_in_output(self):
        bar = _priority_bar(62.5)
        assert "62.5" in bar

    def test_zero_score(self):
        bar = _priority_bar(0.0)
        assert bar  # should not raise or return empty


# ── write_tsv ─────────────────────────────────────────────────────────────────

class TestWriteTsv:
    def test_creates_file(self, tmp_path, sample_collection):
        out = tmp_path / "variants.tsv"
        write_tsv(sample_collection, out)
        assert out.exists()

    def test_header_matches_columns(self, tmp_path, sample_collection):
        out = tmp_path / "variants.tsv"
        write_tsv(sample_collection, out)
        with open(out) as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)
        assert header == TSV_COLUMNS

    def test_row_count(self, tmp_path, sample_collection):
        out = tmp_path / "variants.tsv"
        write_tsv(sample_collection, out)
        with open(out) as f:
            rows = list(csv.reader(f, delimiter="\t"))
        assert len(rows) == len(sample_collection) + 1  # +1 for header

    def test_rank_column_starts_at_one(self, tmp_path, sample_collection):
        out = tmp_path / "variants.tsv"
        write_tsv(sample_collection, out)
        with open(out) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert rows[0]["rank"] == "1"

    def test_none_values_written_as_dot(self, tmp_path):
        v = Variant(chromosome="1", position=100, ref="A", alt="T")
        col = VariantCollection()
        col.add(v)
        out = tmp_path / "variants.tsv"
        write_tsv(col, out)
        with open(out) as f:
            reader = csv.DictReader(f, delimiter="\t")
            row = next(reader)
        assert row["gene"] == "."
        assert row["gnomad_af"] == "."

    def test_accepts_list_input(self, tmp_path, sample_variant):
        out = tmp_path / "list.tsv"
        write_tsv([sample_variant], out)
        assert out.exists()

    def test_creates_parent_directory(self, tmp_path, sample_collection):
        out = tmp_path / "subdir" / "variants.tsv"
        write_tsv(sample_collection, out)
        assert out.exists()


# ── write_html_report ─────────────────────────────────────────────────────────

class TestWriteHtmlReport:
    def test_creates_file(self, tmp_path, sample_collection):
        out = tmp_path / "report.html"
        write_html_report(list(sample_collection), out)
        assert out.exists()

    def test_is_valid_html(self, tmp_path, sample_collection):
        out = tmp_path / "report.html"
        write_html_report(list(sample_collection), out)
        content = out.read_text()
        assert "<!DOCTYPE html>" in content or "<html" in content

    def test_contains_gene_name(self, tmp_path, sample_variant):
        out = tmp_path / "report.html"
        write_html_report([sample_variant], out)
        content = out.read_text()
        assert "BRCA1" in content

    def test_contains_clinvar_badge(self, tmp_path, sample_variant):
        out = tmp_path / "report.html"
        write_html_report([sample_variant], out)
        content = out.read_text()
        assert "#c0392b" in content  # pathogenic red

    def test_contains_omim_link(self, tmp_path, sample_variant):
        out = tmp_path / "report.html"
        write_html_report([sample_variant], out)
        content = out.read_text()
        assert "omim.org" in content

    def test_contains_clinvar_link(self, tmp_path, sample_variant):
        out = tmp_path / "report.html"
        write_html_report([sample_variant], out)
        content = out.read_text()
        assert "ncbi.nlm.nih.gov/clinvar" in content
        assert "12345" in content

    def test_filter_summary_shown(self, tmp_path, sample_variant):
        out = tmp_path / "report.html"
        write_html_report(
            [sample_variant], out,
            filters={"max_af": 0.01, "min_cadd": 15}
        )
        content = out.read_text()
        assert "0.01" in content or "1.0%" in content

    def test_input_filename_in_header(self, tmp_path, sample_variant):
        out = tmp_path / "report.html"
        write_html_report(
            [sample_variant], out,
            input_path=Path("patient.vcf")
        )
        content = out.read_text()
        assert "patient.vcf" in content

    def test_empty_variant_list(self, tmp_path):
        out = tmp_path / "empty.html"
        write_html_report([], out)
        assert out.exists()

    def test_no_external_dependencies(self, tmp_path, sample_collection):
        """HTML report must work offline — no CDN links."""
        out = tmp_path / "report.html"
        write_html_report(list(sample_collection), out)
        content = out.read_text()
        assert "cdn." not in content
        assert "googleapis.com" not in content

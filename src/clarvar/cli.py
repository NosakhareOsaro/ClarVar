"""
Command-line interface for ClarVar.

Your task: wire together the VCFParser, VariantAnnotator, VariantPrioritizer,
and reporter functions into three click commands.

Commands to implement
---------------------

clarvar annotate
    Parses a VCF file and annotates variants using the VEP REST API.
    Writes an annotated VCF to the output directory.

    Options:
        -i / --input    VCF input file (required)
        -o / --output   Output directory (required, created if absent)
        -a / --assembly GRCh38 or GRCh37 (default: GRCh38)
        -v / --verbose  Enable debug logging

clarvar prioritize
    Reads a VCF, scores and ranks variants, writes outputs.
    Writes: prioritised.vcf, prioritised_variants.tsv, report.json, report.html

    Options: same as annotate, plus --no-html and --top-n N

clarvar pipeline
    Runs annotation + prioritisation in one command.
    This is the recommended entry point for most users.
    Writes: prioritised.vcf, prioritised_variants.tsv, report.json, report.html

    Options: same as prioritize

Output files
------------
All commands write into the directory given by -o:

    <output>/annotated.vcf              (annotate only)
    <output>/prioritised.vcf            (prioritize / pipeline)
    <output>/prioritised_variants.tsv   (prioritize / pipeline)
    <output>/report.json                (prioritize / pipeline)
    <output>/report.html                (prioritize / pipeline, unless --no-html)

Implementation notes
--------------------
- Use click.group() for the root CLI and @cli.command() for each subcommand.
- Create the output directory with Path(output).mkdir(parents=True, exist_ok=True).
- Set logging level to DEBUG when --verbose is passed.
- Use click.secho(..., fg="green") for success messages.
- Use click.secho(..., fg="cyan") for step progress inside pipeline.

Example usage (after pip install -e .):

    clarvar pipeline -i data/sample/sample_variants.vcf -o results/
    clarvar pipeline -i patient.vcf -o results/ -a GRCh37 --top-n 100
    clarvar annotate -i raw.vcf -o results/ --verbose
"""

import logging
import json
from pathlib import Path

import click

from .variant import VariantCollection
from .annotator import VariantAnnotator
from .prioritizer import VariantPrioritizer
from .reporter import write_tsv, write_html_report
from .vcf_parser import VCFParser

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
)
logger = logging.getLogger(__name__)


# ── Root command group ────────────────────────────────────────────────────────

@click.group()
@click.version_option(package_name="clarvar")
def cli():
    """ClarVar – Open-Source Variant Annotation & Prioritisation Pipeline."""


# ── annotate ──────────────────────────────────────────────────────────────────

@cli.command()
@click.option("-i", "--input",    required=True, metavar="VCF",
              help="Input VCF file (.vcf or .vcf.gz).")
@click.option("-o", "--output",   required=True, metavar="DIR",
              help="Output directory (created if needed).")
@click.option("-a", "--assembly", default="GRCh38", show_default=True,
              type=click.Choice(["GRCh37", "GRCh38"]),
              help="Genome assembly.")
@click.option("-v", "--verbose",  is_flag=True,
              help="Enable debug logging.")
def annotate(input, output, assembly, verbose):
    """Annotate variants with VEP consequences, gnomAD frequencies, and ClinVar."""
    raise NotImplementedError(
        "TODO: implement annotate command\n\n"
        "Steps:\n"
        "  1. Set logging level to DEBUG if verbose\n"
        "  2. Create output directory\n"
        "  3. Parse the input VCF with VCFParser\n"
        "  4. Annotate with VariantAnnotator(assembly=assembly)\n"
        "  5. Write annotated VCF to <output>/annotated.vcf\n"
        "  6. Print a success message"
    )


# ── prioritize ────────────────────────────────────────────────────────────────

@cli.command()
@click.option("-i", "--input",    required=True, metavar="VCF",
              help="Annotated VCF file.")
@click.option("-o", "--output",   required=True, metavar="DIR",
              help="Output directory (created if needed).")
@click.option("--no-html",  is_flag=True,
              help="Skip HTML report generation.")
@click.option("--top-n",    default=50, show_default=True, metavar="N",
              help="Number of variants shown in the HTML report.")
@click.option("-v", "--verbose",  is_flag=True,
              help="Enable debug logging.")
def prioritize(input, output, no_html, top_n, verbose):
    """Prioritize pre-annotated variants and write TSV, JSON, and HTML reports."""
    raise NotImplementedError(
        "TODO: implement prioritize command\n\n"
        "Steps:\n"
        "  1. Set logging level to DEBUG if verbose\n"
        "  2. Create output directory\n"
        "  3. Parse the input VCF with VCFParser\n"
        "  4. Score and sort with VariantPrioritizer\n"
        "  5. Write prioritised.vcf, prioritised_variants.tsv\n"
        "  6. Write report.json via prioritizer.generate_report()\n"
        "  7. Write report.html via write_html_report() (unless --no-html)\n"
        "  8. Print success messages for each output"
    )


# ── pipeline ──────────────────────────────────────────────────────────────────

@cli.command()
@click.option("-i", "--input",    required=True, metavar="VCF",
              help="Input VCF file (.vcf or .vcf.gz).")
@click.option("-o", "--output",   required=True, metavar="DIR",
              help="Output directory (created if needed).")
@click.option("-a", "--assembly", default="GRCh38", show_default=True,
              type=click.Choice(["GRCh37", "GRCh38"]),
              help="Genome assembly.")
@click.option("--no-html",  is_flag=True,
              help="Skip HTML report generation.")
@click.option("--top-n",    default=50, show_default=True, metavar="N",
              help="Number of variants shown in the HTML report.")
@click.option("-v", "--verbose",  is_flag=True,
              help="Enable debug logging.")
def pipeline(input, output, assembly, no_html, top_n, verbose):
    """Run the full annotation + prioritisation pipeline (recommended)."""
    raise NotImplementedError(
        "TODO: implement pipeline command\n\n"
        "Steps:\n"
        "  1. Set logging level to DEBUG if verbose\n"
        "  2. Create output directory\n"
        "  3. Print a banner header\n"
        "  4. [1/3] Parse VCF\n"
        "  5. [2/3] Annotate with VariantAnnotator(assembly=assembly)\n"
        "  6. [3/3] Score and sort with VariantPrioritizer\n"
        "  7. Write all four output files\n"
        "  8. Print a completion banner\n\n"
        "Hint: you can call the annotate and prioritize logic directly\n"
        "rather than invoking the CLI commands — they share the same steps."
    )


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    cli()


if __name__ == "__main__":
    main()

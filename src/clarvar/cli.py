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
from .hpo_filter import filter_collection_by_hpo

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
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output_dir = Path(output)
    output_dir.mkdir(parents=True, exist_ok=True)

    click.secho("Parsing VCF file...", fg="cyan")

    parser = VCFParser()
    collection = parser.parse_file(input)

    click.secho(f"Annotating {len(collection)} variants using {assembly}...",fg="cyan",)

    annotator = VariantAnnotator(assembly=assembly)

    annotated = annotator.annotate_collection(collection,verbose=verbose,)

    output_file = output_dir / "annotated.vcf"

    parser.write_vcf(annotated,output_file,)

    click.secho(f"Successfully wrote annotated variants to {output_file}",fg="green",)


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
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    out_dir = Path(output)
    out_dir.mkdir(parents=True, exist_ok=True)
    parser = VCFParser()
    variants = parser.parse_file(Path(input))
    prioritizer = VariantPrioritizer()
    prioritized = prioritizer.prioritize_collection(variants)
    parser.write_vcf(prioritized, out_dir / "prioritised.vcf")
    write_tsv(prioritized, out_dir / "prioritised_variants.tsv")
    (out_dir / "report.json").write_text(json.dumps(prioritizer.generate_report(prioritized), indent=2))
    if not no_html:
        write_html_report(list(prioritized.variants)[:top_n], out_dir / "report.html", input_path=Path(input))
    click.secho("✓ Done.", fg="green")


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
@click.option("--hpo", multiple=True,
              help="Filter output to genes associated with one or more HPO terms.")
@click.option("-v", "--verbose",  is_flag=True,
              help="Enable debug logging.")
def pipeline(input, output, assembly, no_html, top_n, hpo, verbose):
    """Run the full annotation + prioritisation pipeline (recommended)."""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    out_dir = Path(output)
    out_dir.mkdir(parents=True, exist_ok=True)
    click.echo("=" * 55)
    click.echo("  ClarVar: Variant Annotation & Prioritisation")
    click.echo("=" * 55)
    click.secho(f"\n[1/3] Parsing   {input}", fg="cyan")
    parser = VCFParser()
    variants = parser.parse_file(Path(input))
    click.secho(f"      ✓ {len(variants)} variants parsed", fg="cyan")
    click.secho("\n[2/3] Annotating ...", fg="cyan")
    annotator = VariantAnnotator(assembly=assembly)
    annotated = annotator.annotate_collection(variants, verbose=verbose)
    click.secho(f"      ✓ {len(annotated)} variants annotated", fg="cyan")
    click.secho("\n[3/3] Prioritising ...", fg="cyan")
    prioritizer = VariantPrioritizer()
    prioritized = prioritizer.prioritize_collection(annotated)
    if hpo:
        prioritized = filter_collection_by_hpo(prioritized, list(hpo))
    click.secho(f"      ✓ {len(prioritized)} variants scored", fg="cyan")
    parser.write_vcf(prioritized, out_dir / "prioritised.vcf")
    write_tsv(prioritized, out_dir / "prioritised_variants.tsv")
    (out_dir / "report.json").write_text(json.dumps(prioritizer.generate_report(prioritized), indent=2))
    if not no_html:
        write_html_report(list(prioritized.variants)[:top_n], out_dir / "report.html", input_path=Path(input))
    click.secho("\n  Pipeline complete!", fg="green", bold=True)
    click.echo("=" * 55)


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    cli()


if __name__ == "__main__":
    main()

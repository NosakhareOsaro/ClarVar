"""
Reporter module — writes prioritised variants to TSV and HTML output files.

Your task: implement two output formats.

TSV output
----------
A tab-separated file with one variant per row and the following columns
(in this order):

    rank, chrom, pos, ref, alt, gene, consequence, hgvsc, hgvsp,
    gnomad_af, cadd_phred, clinvar_sig, clinvar_id, sift, polyphen,
    priority_score

Use "." for any field that is None.
Format floats: gnomad_af to 6 decimal places, all others to 2.

HTML report
-----------
A self-contained HTML file (no external JS/CSS dependencies) that can be
opened in any browser. It should include:

1. A header with the input filename and generation timestamp.
2. A filter summary bar showing which gnomAD AF / CADD / consequence
   filters were applied (pass these in via the `filters` dict).
3. A sortable/searchable HTML table with columns matching the TSV, plus:
   - Colour-coded ClinVar badges:
       pathogenic          → red    (#c0392b)
       likely_pathogenic   → orange (#e67e22)
       uncertain_sig       → purple (#8e44ad)
       benign/likely_benign → green (#27ae60)
   - A visual priority score bar (a coloured div scaled to the score)
   - Clickable links:
       gene name  → https://www.omim.org/search?search=<GENE>
       clinvar_id → https://www.ncbi.nlm.nih.gov/clinvar/variation/<ID>/
4. A footer crediting the data sources.

Tip: build the HTML as an f-string. All user-supplied strings should be
passed through html.escape() before inserting into HTML to prevent injection.

Example usage
-------------
>>> from clarvar.reporter import write_tsv, write_html_report
>>> from pathlib import Path
>>>
>>> write_tsv(prioritized, Path("results/variants.tsv"))
>>> write_html_report(
...     list(prioritized)[:50],
...     Path("results/report.html"),
...     input_path=Path("patient.vcf"),
...     filters={"max_af": 0.01, "min_cadd": 15},
... )
"""

import csv
import html as html_module
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Union

from clarvar.variant import Variant, VariantCollection

# Column names for the TSV output — do not change this order
TSV_COLUMNS = [
    "rank", "chrom", "pos", "ref", "alt",
    "gene", "consequence", "hgvsc", "hgvsp",
    "gnomad_af", "cadd_phred",
    "clinvar_sig", "clinvar_id",
    "sift", "polyphen",
    "priority_score",
]


def write_tsv(
    variants: Union[List[Variant], VariantCollection],
    path: Path,
) -> None:
    """
    Write prioritised variants to a tab-separated file.

    Columns are defined in TSV_COLUMNS. Variants are written in the order
    they appear in the input (assumed already sorted by priority_score).
    Missing values (None) should be written as ".".

    Parameters
    ----------
    variants : list of Variant or VariantCollection
    path : Path
        Output file path. Parent directories are created if needed.

    Example
    -------
    >>> write_tsv(prioritized, Path("results/variants.tsv"))
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("\t".join(TSV_COLUMNS) + "\n")
        for rank, variant in enumerate(variants, start=1):
            row = []
            for col in TSV_COLUMNS:
                precision = 6 if col == "gnomad_af" else 2
                if col == "rank":
                    value = rank
                else:
                    attr = COLUMN_ATTR_MAP.get(col, col)  # fall back to col itself if no mapping needed
                    value = getattr(variant, attr)
                row.append(_fmt(value, precision))
            f.write("\t".join(row) + "\n")


def write_html_report(
    variants: Union[List[Variant], VariantCollection],
    path: Path,
    input_path: Optional[Path] = None,
    filters: Optional[Dict] = None,
    title: str = "ClarVar — Variant Prioritisation Report",
) -> None:
    """
    Generate a self-contained HTML report for a list of prioritised variants.

    The report should be a single .html file that works without internet
    access (no CDN links). All styling should be inline or in a <style> block.

    Parameters
    ----------
    variants : list of Variant or VariantCollection
        The variants to display (typically the top N from a prioritised list).
    path : Path
        Output file path.
    input_path : Path or None
        Original VCF input path — used in the report header.
    filters : dict or None
        Filter settings used during prioritisation, e.g.:
        {"max_af": 0.01, "min_cadd": 15, "consequences": {"missense_variant"}}
        Displayed in the filter summary bar.
    title : str
        HTML page title and report heading.

    Example
    -------
    >>> write_html_report(
    ...     list(prioritized)[:50],
    ...     Path("results/report.html"),
    ...     input_path=Path("patient.vcf"),
    ...     filters={"max_af": 0.01, "min_cadd": 15},
    ... )
    """

    # report title
    title_html = f"<h1>{html_module.escape(title)}</h1>"

    # add input_path as optional header subtitle
    input_html = ""
    if input_path:
        input_html = f'<div style="color:#7f8c8d;font-size:0.9em;margin-bottom:10px;">Input: {html_module.escape(str(input_path))}</div>'

    # generation timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    timestamp_html = f'<div style="color:#7f8c8d;font-size:0.9em;margin-bottom:20px;">Generated on {timestamp}</div>'

    # report filter summary bar
    if filters:
        filter_items = []
        if "max_af" in filters:
            filter_items.append(f"gnomAD AF &le; {filters['max_af']}")
        if "min_cadd" in filters:
            filter_items.append(f"CADD &ge; {filters['min_cadd']}")
        if "consequences" in filters:
            cons = ", ".join(filters["consequences"])
            filter_items.append(f"Consequences: {cons}")
        filter_summary = " | ".join(filter_items)
        filter_html = f'<div style="background:#ecf0f1;padding:10px;border-radius:5px;margin-bottom:20px;">{filter_summary}</div>'
    else:
        filter_html = ""

    # color coded table rows with ClinVar badges and priority bars
    table_rows = []
    for var in variants:
        clinvar_badge = _clinvar_badge(var.clinvar_significance) 
        priority_bar = _priority_bar(var.priority_score)
        row_html = f"""
        <tr>
            <td>{var.chromosome}:{var.position} {var.ref}>{var.alt}</td>
            <td><a href="https://www.omim.org/search?search={html_module.escape(var.gene or '')}">{html_module.escape(var.gene or '.')}</a></td>
            <td>{html_module.escape(str(var.consequence) if var.consequence else ".")}</td>
            <td>{html_module.escape(var.hgvsc or ".")}</td>
            <td>{html_module.escape(var.hgvsp or ".")}</td>
            <td>{_fmt(var.allele_frequency, 6)}</td>
            <td>{_fmt(var.cadd_phred)}</td>
            <td>{clinvar_badge}</td>
            <td>{f'<a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/{html_module.escape(var.clinvar_id)}/">{html_module.escape(var.clinvar_id)}</a>' if var.clinvar_id else "."}</td>
            <td>{_fmt(var.sift)}</td>
            <td>{_fmt(var.polyphen)}</td>
            <td>{priority_bar}</td>
        </tr>
        """
        table_rows.append(row_html)
    
    table_html = f"""
    <table style="width:100%;border-collapse:collapse;">
        <thead>
            <tr>
                <th>Variant</th>
                <th>Gene</th>
                <th>Consequence</th>
                <th>HGVSc</th>
                <th>HGVSp</th>
                <th>gnomAD AF</th>
                <th>CADD</th>
                <th>ClinVar Sig</th>
                <th>ClinVar ID</th>
                <th>SIFT</th>
                <th>PolyPhen</th>
                <th>Priority Score</th>
            </tr>
        </thead>
        <tbody>
            {''.join(table_rows)}
        </tbody>
    </table>
    """

    # generate final HTML
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>{html_module.escape(title)}</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            table, th, td {{ border: 1px solid #bdc3c7; }}
            th, td {{ padding: 8px; text-align: left; }}
            th {{ background-color: #ecf0f1; }}
            tr:nth-child(even) {{ background-color: #f9f9f9; }}
        </style>
    </head>
    <body>
        {title_html}
        {input_html}
        {timestamp_html}
        {filter_html}
        {table_html}
    </body>
    </html>
    """

    # write to file (added date-based suffix if file already exists to avoid overwriting)
    if path.exists():
        timestamp_suffix = datetime.now().strftime("%Y%m%d%H%M%S")
        path = path.with_name(f"{path.stem}_{timestamp_suffix}{path.suffix}")
        print(f"Warning: report already exists. Writing to {path} instead.")

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w",  encoding="utf-8") as f:
        f.write(html_content)

    return path


# ── Internal helpers — implement these to keep write_html_report clean ────────

def _fmt(value, precision: int = 2) -> str:
    """
    Format a value for table display.

    Return "." for None. Format floats to `precision` decimal places.
    Return str(value) for everything else.

    Parameters
    ----------
    value : any
    precision : int

    Returns
    -------
    str
    """
    if value is None or value == "None":
        return "."
    elif isinstance(value, float):
        return f"{value:.{precision}f}"
    else:
        return str(value)


def _clinvar_badge(significance: Optional[str]) -> str:
    """
    Return an HTML <span> badge styled by ClinVar significance tier.

    Return an empty string if significance is None.

    Colour mapping:
        pathogenic            → #c0392b (red)
        likely_pathogenic     → #e67e22 (orange)
        uncertain_significance → #8e44ad (purple)
        likely_benign / benign → #27ae60 (green)
        anything else          → #7f8c8d (grey)

    The badge should be a white-text pill:
    <span style="background:{colour};color:#fff;padding:2px 6px;
                 border-radius:3px;font-size:0.8em;">{significance}</span>

    Parameters
    ----------
    significance : str or None

    Returns
    -------
    str
        HTML string for the badge, or "" if significance is None.
    """
    colours = {
        "pathogenic": "#c0392b",
        "likely_pathogenic": "#e67e22",
        "uncertain_significance": "#8e44ad",
        "likely_benign": "#27ae60",
        "benign": "#27ae60",
    }
    
    if not significance:
        return ""

    sig = significance.lower()
    colour = colours.get(sig, "#7f8c8d")
    label = html_module.escape(sig)
    
    return f'<span style="background:{colour}; color:#fff; padding:2px 6px; border-radius:3px;">{label}</span>'


def _priority_bar(score: float, max_score: float = 100) -> str:
    """
    Return an HTML snippet showing a coloured progress bar for a priority score.

    The bar should be proportional to score/max_score (capped at 100%).
    Suggested colour thresholds:
        score >= 60 → red    (#e74c3c)
        score >= 35 → orange (#e67e22)
        else        → blue   (#3498db)

    The bar should sit inside a fixed-width grey background div, with the
    score value displayed as text to the right.

    Parameters
    ----------
    score : float
    max_score : float

    Returns
    -------
    str
        HTML string.
    """

    if score is None:
        score = 0.0

    # Ensure numeric and round for display
    try:
        score = float(score)
    except Exception:
        score = 0.0

    score_display = round(score, 1)

    # determine width as percentage of max_score, capped at 100%
    if max_score <= 0:
        pct = 0
    else:
        pct = max(0.0, min(100.0, (score / max_score) * 100.0))

    # colour changes based on severity
    if score >= 60:
        colour = "#e74c3c"
    elif score >= 35:
        colour = "#e67e22"
    else:
        colour = "#3498db"

    # fixed-width background with inner coloured bar and score text to the right
    return (
        f'<div style="display:inline-flex; align-items:center; gap:8px;">'
        f'<div style="background:#eee; border-radius:4px; width:100px; height:16px; overflow:hidden;">'
        f'<div style="background:{colour}; width:{pct}%; height:100%;"></div>'
        f'</div>'
        f'<span style="font-size:0.9em;">{score_display}</span>'
        f'</div>'
    )

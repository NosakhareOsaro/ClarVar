from pathlib import Path
from typing import Dict, Iterable, Set
import csv
import urllib.request

from .variant import VariantCollection

HPO_URL = "https://hpo.jax.org/app/data/annotations/genes_to_phenotype.txt"

def get_cache_path() -> Path:
    """Return the local cache path for the HPO annotation file."""
    cache_dir = Path.home() / ".cache" / "clarvar"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / "genes_to_phenotype.txt"

def download_hpo_annotations(force: bool = False) -> Path:
    """Download and cache the HPO gene-to-phenotype annotation file."""
    cache_path = get_cache_path()
    
    if cache_path.exists() and not force:
        return cache_path

    urllib.request.urlretrieve(HPO_URL, cache_path)
    return cache_path

def build_hpo_gene_lookup(annotation_file: Path | None = None) -> Dict[str, Set[str]]:
    """Build a lookup mapping HPO term IDs to associated gene symbols."""
    if annotation_file is None:
        annotation_file = download_hpo_annotations()

    lookup: Dict[str, Set[str]] = {}

    with annotation_file.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            hpo_id = row.get("hpo_term_id")
            gene_symbol = row.get("entrez_gene_symbol")

            if not hpo_id or not gene_symbol:
                continue

            lookup.setdefault(hpo_id, set()).add(gene_symbol)

    return lookup

def genes_for_hpo_terms(hpo_terms: Iterable[str]) -> Set[str]:
    """Return all genes associated with the provided HPO terms."""
    lookup = build_hpo_gene_lookup()
    genes: Set[str] = set()

    for term in hpo_terms:
        genes.update(lookup.get(term, set()))

    return genes


def filter_collection_by_hpo(
    collection: VariantCollection,
    hpo_terms: Iterable[str],
) -> VariantCollection:
    """Filter variants to genes associated with the provided HPO terms."""
    allowed_genes = genes_for_hpo_terms(hpo_terms)

    filtered_variants = [
        variant
        for variant in collection
        if variant.gene in allowed_genes
    ]

    return VariantCollection(
        variants=filtered_variants,
        name=f"{collection.name}_hpo_filtered",
    )


    


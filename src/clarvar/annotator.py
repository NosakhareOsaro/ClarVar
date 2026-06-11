"""
Variant annotation module.

Your task: enrich each Variant object with real data from external sources.

Data sources to integrate
--------------------------
1. Ensembl VEP REST API  — consequence terms, CADD PHRED, SIFT, PolyPhen,
                           HGVSc/HGVSp notation, canonical transcript
   GRCh38 endpoint:  https://rest.ensembl.org/vep/human/region
   GRCh37 endpoint:  https://grch37.rest.ensembl.org/vep/human/region
   Docs:             https://rest.ensembl.org/documentation/info/vep_region_post

2. gnomAD v4.1         — population allele frequencies
   Returned inside VEP colocated_variants → frequencies

3. ClinVar             — clinical significance classifications
   Returned inside VEP colocated_variants → clin_sig

API notes
---------
- VEP accepts up to 200 variants per POST request (see BATCH_SIZE).
- Variants must be formatted as: "chrom pos end ref/alt 1"
  e.g. SNV:       "1 925952 925952 G/A 1"
       deletion:  "7 117548628 117548630 CTT/C 1"  (end = pos + len(ref) - 1)
       insertion: "2 47702543 47702543 C/CT 1"
- The API returns HTTP 429 when rate-limited; respect the Retry-After header.
- Request both CADD and gnomAD annotations: include "CADD": 1, "gnomAD": 1,
  "ClinVar": 1 in the POST body.

VEP response structure (simplified)
-------------------------------------
[
  {
    "input": "1 925952 925952 G/A 1",
    "transcript_consequences": [
      {
        "gene_symbol": "SAMD11",
        "transcript_id": "ENST00000342066",
        "consequence_terms": ["missense_variant"],
        "hgvsc": "ENST00000342066.3:c.599G>A",
        "hgvsp": "ENSP00000342992.3:p.Arg200Gln",
        "sift_prediction": "deleterious",
        "polyphen_prediction": "probably_damaging",
        "cadd_phred": 24.1
      },
      ...
    ],
    "colocated_variants": [
      {
        "clin_sig": ["likely_benign"],
        "frequencies": {
          "A": { "gnomad": 0.0003, "gnomadg": 0.0002 }
        }
      }
    ]
  },
  ...
]

When multiple transcript_consequences are returned, pick the most severe one
using the VEP consequence severity order (see CONSEQUENCE_SEVERITY below).
"""

import logging
import time
from typing import Any, Dict, List, Optional

import requests

from clarvar.variant import Consequence, Variant, VariantCollection

logger = logging.getLogger(__name__)

# ── Constants — do not change these ──────────────────────────────────────────

VEP_ENDPOINTS: Dict[str, str] = {
    "GRCh38": "https://rest.ensembl.org/vep/human/region",
    "GRCh37": "https://grch37.rest.ensembl.org/vep/human/region",
}

BATCH_SIZE = 200        # maximum variants per VEP POST request
RATE_LIMIT_PAUSE = 0.5  # seconds to wait between batches
MAX_RETRIES = 3         # retry attempts on failure

# Lower rank = more severe consequence
CONSEQUENCE_SEVERITY: Dict[str, int] = {
    "transcript_ablation": 1,
    "splice_acceptor_variant": 2,
    "splice_donor_variant": 3,
    "stop_gained": 4,
    "frameshift_variant": 5,
    "stop_lost": 6,
    "start_lost": 7,
    "transcript_amplification": 8,
    "inframe_insertion": 9,
    "inframe_deletion": 10,
    "disruptive_inframe_deletion": 11,
    "missense_variant": 12,
    "splice_region_variant": 14,
    "synonymous_variant": 17,
    "upstream_gene_variant": 26,
    "downstream_gene_variant": 27,
    "intron_variant": 23,
    "intergenic_variant": 36,
}


# ── Helper functions — implement these first ──────────────────────────────────

def _to_vep_region(variant: Variant) -> str:
    """
    Convert a Variant to the Ensembl VEP region string format.

    Format: "{chrom} {pos} {end} {ref}/{alt} 1"
    where end = pos + len(ref) - 1

    Examples
    --------
    SNV   chr1:100 A>T   → "1 100 100 A/T 1"
    Del   chr7:100 ATG>A → "7 100 102 ATG/A 1"
    Ins   chr2:100 A>ATG → "2 100 100 A/ATG 1"

    Parameters
    ----------
    variant : Variant

    Returns
    -------
    str
        VEP region notation string.
    """
    
    end = variant.position + len(variant.ref) - 1
    return f"{variant.chromosome} {variant.position} {end} {variant.ref}/{variant.alt} 1"

def _pick_most_severe_transcript(
    transcripts: List[Dict[str, Any]]
) -> Optional[Dict[str, Any]]:
    """
    Return the transcript consequence dict with the most severe consequence term.

    Use CONSEQUENCE_SEVERITY to rank terms. If a transcript has multiple
    consequence_terms, use the most severe one for ranking. Return the
    transcript dict itself (not just the term), so the caller can access
    all fields (gene_symbol, hgvsc, cadd_phred, etc.).

    Parameters
    ----------
    transcripts : list of dicts
        The "transcript_consequences" list from a VEP response hit.

    Returns
    -------
    dict or None
        The most severe transcript dict, or None if the list is empty.
    """
    if not transcripts:
        return None

    most_severe = None
    lowest_rank = float('inf')

    for transcript in transcripts:
        for term in transcript.get("consequence_terms", []):
            rank = CONSEQUENCE_SEVERITY.get(term, float('inf'))
            if rank < lowest_rank:
                lowest_rank = rank
                most_severe = transcript

    return most_severe


def _extract_gnomad_af(vep_hit: Dict[str, Any]) -> Optional[float]:
    """
    Extract the gnomAD allele frequency from a VEP response hit.

    gnomAD frequencies live in:
      vep_hit["colocated_variants"][i]["frequencies"][allele]["gnomad"]
    or under keys "gnomadg" (genomes) or "gnomade" (exomes).

    Return the first gnomAD frequency found, or None if absent.

    Parameters
    ----------
    vep_hit : dict
        A single element from the VEP JSON response list.

    Returns
    -------
    float or None
    """
    colocated = vep_hit.get("colocated_variants", [])

    for variant in colocated:
        frequencies = variant.get("frequencies", {})

        for allele_data in frequencies.values():
            for key in ("gnomad", "gnomadg", "gnomade"):
                if key in allele_data:
                    return allele_data[key]

    return None

def _extract_clinvar(
    vep_hit: Dict[str, Any],
) -> tuple[Optional[str], Optional[str]]:
    """
    Extract ClinVar clinical significance and variant ID from a VEP hit.

    ClinVar data lives in colocated_variants:
      vep_hit["colocated_variants"][i]["clin_sig"]       → list of strings
      vep_hit["colocated_variants"][i]["var_synonyms"]["ClinVar"] → list

    Skip entries where "somatic" is True.

    Parameters
    ----------
    vep_hit : dict

    Returns
    -------
    (significance, clinvar_id) : tuple of (str or None, str or None)
        significance — comma-joined string if multiple values
        clinvar_id   — first ClinVar ID, or None
    """
    for colocated in vep_hit.get("colocated_variants", []):
        # Skip somatic (e.g. tumour-acquired) records — ClarVar is germline-only.
        if colocated.get("somatic"):
            continue

        clin_sig = colocated.get("clin_sig")
        if clin_sig:
            significance = ", ".join(clin_sig) if isinstance(clin_sig, list) else str(clin_sig)
            clinvar_ids = colocated.get("var_synonyms", {}).get("ClinVar")
            clinvar_id = clinvar_ids[0] if clinvar_ids else None
            return significance, clinvar_id

    return None, None


def _apply_vep_result(variant: Variant, vep_hit: Dict[str, Any]) -> None:
    """
    Write VEP annotation fields back onto a Variant object in-place.

    Fields to populate from the most severe transcript:
      variant.gene        ← gene_symbol (or gene_id as fallback)
      variant.transcript  ← transcript_id
      variant.consequence ← Consequence.from_vep_string(consequence_terms[0])
      variant.hgvsc       ← hgvsc
      variant.hgvsp       ← hgvsp
      variant.sift        ← sift_prediction
      variant.polyphen    ← polyphen_prediction
      variant.cadd_phred  ← cadd_phred (cast to float)

    Fields to populate from colocated_variants:
      variant.allele_frequency    ← from _extract_gnomad_af()
      variant.clinvar_significance,
      variant.clinvar_id          ← from _extract_clinvar()

    Parameters
    ----------
    variant : Variant
        Modified in-place.
    vep_hit : dict
        A single element from the VEP JSON response list.
    """
    transcript=_pick_most_severe_transcript(vep_hit)
    if transcript is None:
        return
    transcripts = vep_hit.get("transcript_consequences", [])
    transcript = _pick_most_severe_transcript(transcripts)
    
    variant.gene=(transcript.get("gene_symbol") or transcript.get("gene_id"))
             
    variant.transcript=transcript.get("transcript_id")

    terms = transcript.get("consequence_terms", [])
    if terms:
        variant.consequence = Consequence.from_vep_string(terms[0])

    variant.hgvsc = transcript.get("hgvsc")
    variant.hgvsp = transcript.get("hgvsp")

    variant.sift = transcript.get("sift_prediction")
    variant.polyphen = transcript.get("polyphen_prediction")
   
    cadd= transcript.get("cadd_phred")
    if cadd is not None:
        variant.cadd_phred=float(cadd)
    
    variant.allele_frequency = _extract_gnomad_af(vep_hit)

    sig, cid = _extract_clinvar(vep_hit)
    variant.clinvar_significance = sig
    variant.clinvar_id = cid


def _post_vep_batch(
    variants: List[Variant],
    endpoint: str,
) -> List[Dict[str, Any]]:
    """
    POST a batch of variants to the Ensembl VEP REST endpoint.

    Build a JSON payload like:
    {
        "variants": ["1 100 100 A/T 1", ...],
        "CADD": 1,
        "gnomAD": 1,
        "ClinVar": 1,
        "canonical": 1,
        "pick": 1
    }

    Handle:
    - HTTP 200 → return parsed JSON list
    - HTTP 429 → wait Retry-After seconds, then retry
    - Other errors / timeouts → retry up to MAX_RETRIES times,
      then log a warning and return []

    Parameters
    ----------
    variants : list of Variant
    endpoint : str
        The full VEP REST URL.

    Returns
    -------
    list of dict
        Parsed JSON response, or [] on failure.
    """
    payload = {
        "variants": [_to_vep_region(v) for v in variants],
        "CADD": 1,
        "gnomAD": 1,
        "ClinVar": 1,
        "canonical": 1,
        "pick": 1,
    }
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    for attempt in range(MAX_RETRIES):
        try:
            response = requests.post(
                endpoint, json=payload, headers=headers, timeout=30
            )

            if response.status_code == 200:
                return response.json()

            if response.status_code == 429:
                wait = int(response.headers.get("Retry-After", 1))
                logger.warning("VEP rate-limited; waiting %ss before retry", wait)
                time.sleep(wait)
                continue

            logger.warning("VEP returned HTTP %s", response.status_code)

        except requests.RequestException as err:
            logger.warning("VEP request failed: %s", err)

        time.sleep(2 ** attempt)

    logger.warning("VEP batch failed after %s attempts; returning []", MAX_RETRIES)
    return []


# ── Public classes ────────────────────────────────────────────────────────────

class VariantAnnotator:
    """
    Annotates variants using the Ensembl VEP REST API.

    Sends variants in batches of BATCH_SIZE, respects rate limits,
    and writes annotation results directly onto each Variant object.

    Parameters
    ----------
    assembly : str
        "GRCh38" (default) or "GRCh37".
    timeout : int
        HTTP timeout in seconds.
    cache : bool
        Reserved for future caching support.

    Example
    -------
    >>> annotator = VariantAnnotator(assembly="GRCh38")
    >>> annotated = annotator.annotate_collection(variants, verbose=True)
    """

    def __init__(
        self,
        assembly: str = "GRCh38",
        timeout: int = 30,
        cache: bool = True,
    ):
        self.assembly = assembly
        self.timeout = timeout
        self.cache = cache
        self._endpoint = VEP_ENDPOINTS.get(assembly, VEP_ENDPOINTS["GRCh38"])
        self._response_cache: Dict[str, Any] = {}

    def annotate_variant(self, variant: Variant) -> Variant:
        """
        Annotate a single Variant via the VEP REST API.

        Convenience wrapper around annotate_collection() for single variants.
        Prefer annotate_collection() when you have many variants.

        Parameters
        ----------
        variant : Variant

        Returns
        -------
        Variant
            The same object, with annotation fields populated in-place.
        """
        raise NotImplementedError("TODO: implement annotate_variant")

    def annotate_collection(
        self,
        collection: VariantCollection,
        verbose: bool = False,
    ) -> VariantCollection:
        """
        Annotate all variants in a VariantCollection via the VEP REST API.

        Sends variants in batches of BATCH_SIZE. After each batch, waits
        RATE_LIMIT_PAUSE seconds to respect Ensembl's rate limit.

        Implementation steps:
        1. Build a dict mapping VEP region string → Variant object
        2. Loop over variants in slices of BATCH_SIZE
        3. Call _post_vep_batch() for each slice
        4. For each hit in the response, look up the matching Variant by
           hit["input"] and call _apply_vep_result()
        5. Sleep RATE_LIMIT_PAUSE between batches
        6. Return a new VariantCollection with the annotated variants

        Parameters
        ----------
        collection : VariantCollection
        verbose : bool
            If True, print progress to stderr.

        Returns
        -------
        VariantCollection
            New collection with annotation fields populated on every Variant.
        """
        raise NotImplementedError("TODO: implement annotate_collection")

    def clear_cache(self) -> None:
        """Clear the internal response cache."""
        self._response_cache.clear()
        logger.info("Annotation cache cleared")


class LocalAnnotator:
    """
    Lightweight local annotation using pre-downloaded database files.

    Useful for offline annotation or HPC environments without internet access.
    Populate gnomad_db and clinvar_db with variant_key → data mappings,
    where variant_key = "{chrom}:{pos}_{ref}_{alt}".

    Example
    -------
    >>> la = LocalAnnotator(db_path="/data/annotation_db")
    >>> la.load_databases()
    >>> annotated = la.annotate_collection(variants)
    """

    def __init__(self, db_path: str = "./data/annotation_db"):
        self.db_path = db_path
        self.gnomad_db: Dict[str, Dict[str, Any]] = {}
        self.clinvar_db: Dict[str, Dict[str, Any]] = {}

    def load_databases(self) -> None:
        """
        Load annotation databases from local files.

        Implement this to read pre-downloaded gnomAD and ClinVar TSV/VCF
        files into self.gnomad_db and self.clinvar_db.
        """
        logger.info(f"Loading databases from {self.db_path}")
        # TODO: implement database loading from local files

    def annotate_variant(self, variant: Variant) -> Variant:
        """
        Annotate a single variant using the local databases.

        Looks up the variant key "{chrom}:{pos}_{ref}_{alt}" in gnomad_db
        and clinvar_db and populates allele_frequency and
        clinvar_significance if found.

        Parameters
        ----------
        variant : Variant

        Returns
        -------
        Variant
            Same object, modified in-place.
        """
        raise NotImplementedError("TODO: implement LocalAnnotator.annotate_variant")

    def annotate_collection(self, collection: VariantCollection) -> VariantCollection:
        """
        Annotate all variants in a collection using local databases.

        Parameters
        ----------
        collection : VariantCollection

        Returns
        -------
        VariantCollection
        """
        raise NotImplementedError("TODO: implement LocalAnnotator.annotate_collection")

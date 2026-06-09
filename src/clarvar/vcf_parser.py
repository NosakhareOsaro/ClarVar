"""
VCF (Variant Call Format) file parser.

Reads standard VCF format files (plain .vcf or gzip-compressed .vcf.gz),
returning a VariantCollection. Handles multi-allelic sites (takes first ALT),
skips structural variants (<DEL> etc.), and normalises chromosome names
to the no-'chr' prefix convention used by Ensembl VEP.
"""

import gzip
import logging
import re
from pathlib import Path
from typing import List, Optional

from clarvar.variant import Variant, VariantCollection, Consequence

logger = logging.getLogger(__name__)

VCF_VERSION = "VCFv4.2"


def _open_vcf(path: Path):
    """Open a plain or gzip-compressed VCF file for reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _parse_info(info_str: str) -> dict:
    """Parse the INFO column string into a dict.

    Flag keys (no '=') are stored as True.
    """
    if info_str in (".", ""):
        return {}
    result = {}
    for token in info_str.split(";"):
        if "=" in token:
            key, _, value = token.partition("=")
            result[key] = value
        else:
            result[token] = True
    return result


class VCFParser:
    """
    Parser for VCF (Variant Call Format) files.

    VCF specification: https://samtools.github.io/hts-specs/VCFv4.2.pdf
    """

    def __init__(self, validate: bool = True):
        """
        Parameters
        ----------
        validate : bool
            Whether to warn on malformed lines (currently always True).
        """
        self.validate = validate
        self.sample_names: List[str] = []

    def parse_file(self, filepath: str | Path) -> VariantCollection:
        """
        Parse a VCF file and return a VariantCollection.

        Supports plain .vcf and gzip-compressed .vcf.gz.
        Skips structural variants (ALT starts with '<').
        For multi-allelic sites, only the first ALT allele is used.
        Chromosome names are normalised to remove any 'chr' prefix.

        Parameters
        ----------
        filepath : str or Path
            Path to the VCF file.

        Returns
        -------
        VariantCollection
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"VCF file not found: {filepath}")

        logger.info(f"Parsing VCF file: {filepath}")
        collection = VariantCollection(name=filepath.stem)
        count = 0

        with _open_vcf(filepath) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    self._parse_header_line(line)
                    continue
                if line.startswith("#"):
                    continue

                variant = self._parse_variant_line(line)
                if variant:
                    collection.add(variant)
                    count += 1
                    if count % 1000 == 0:
                        logger.debug(f"Parsed {count} variants …")

        logger.info(f"Parsed {count} variants from {filepath}")
        return collection

    def parse_lines(self, lines: List[str]) -> VariantCollection:
        """
        Parse a list of VCF-format strings without file I/O.

        Useful for testing.
        """
        collection = VariantCollection(name="variants")
        for line in lines:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    self._parse_header_line(line)
                continue
            variant = self._parse_variant_line(line.strip())
            if variant:
                collection.add(variant)
        logger.info(f"Parsed {len(collection)} variants from {len(lines)} lines")
        return collection

    def _parse_header_line(self, line: str) -> None:
        """Extract sample names from the #CHROM header line."""
        fields = line.strip().split("\t")
        if len(fields) > 9:
            self.sample_names = fields[9:]
            logger.debug(f"Found {len(self.sample_names)} sample(s): {self.sample_names}")

    def _parse_variant_line(self, line: str) -> Optional[Variant]:
        """
        Parse a single VCF data line into a Variant.

        Returns None for structural variants or malformed lines.
        """
        if not line or line.startswith("#"):
            return None

        try:
            fields = line.split("\t")
            if len(fields) < 8:
                logger.warning(f"Malformed VCF line (too few fields): {line[:80]}")
                return None

            chrom, pos_str, vid, ref, alt_field, qual_str, filt, info_str = fields[:8]

            # Skip structural variants
            first_alt = alt_field.split(",")[0]
            if re.match(r"<[A-Z]+>", first_alt):
                return None

            # Normalise chromosome — strip 'chr' prefix for Ensembl VEP compat
            chrom = chrom.lstrip("chr")

            pos = int(pos_str)

            qual: Optional[float] = None
            if qual_str not in (".", ""):
                try:
                    qual = float(qual_str)
                except ValueError:
                    pass

            # Multi-allelic: take only first ALT
            alt = first_alt

            info = _parse_info(info_str)

            # Optional genotype
            sample_id = "."
            genotype = "."
            if len(fields) >= 10:
                sample_id = self.sample_names[0] if self.sample_names else fields[9]
                if ":" in fields[9]:
                    genotype = fields[9].split(":")[0]
                else:
                    genotype = fields[9]

            return Variant(
                chromosome=chrom,
                position=pos,
                ref=ref,
                alt=alt,
                sample_id=sample_id,
                genotype=genotype,
                quality=qual or 0.0,
                annotation=info,
            )

        except (ValueError, IndexError) as exc:
            logger.warning(f"Failed to parse VCF line ({exc}): {line[:80]}")
            return None

    def write_vcf(self, collection: VariantCollection, filepath: str | Path) -> None:
        """
        Write an annotated VariantCollection back to VCF format.

        Annotation fields (consequence, AF, CADD, ClinVar, priority score)
        are encoded in the INFO column.
        """
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        logger.info(f"Writing {len(collection)} variants to {filepath}")

        with open(filepath, "w", encoding="utf-8") as fh:
            fh.write(f"##fileformat={VCF_VERSION}\n")
            fh.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP consequence">\n')
            fh.write('##INFO=<ID=AF,Number=A,Type=Float,Description="gnomAD allele frequency">\n')
            fh.write('##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD PHRED score">\n')
            fh.write('##INFO=<ID=CLINVAR,Number=1,Type=String,Description="ClinVar significance">\n')
            fh.write('##INFO=<ID=PRIORITY,Number=1,Type=Float,Description="ClarVar priority score">\n')
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

            for v in collection:
                info_parts = []
                if v.consequence and v.consequence != Consequence.UNKNOWN:
                    info_parts.append(f"CSQ={v.consequence.value}")
                if v.allele_frequency is not None:
                    info_parts.append(f"AF={v.allele_frequency:.6f}")
                if v.cadd_phred is not None:
                    info_parts.append(f"CADD={v.cadd_phred:.2f}")
                if v.clinvar_significance:
                    info_parts.append(f"CLINVAR={v.clinvar_significance}")
                if v.priority_score > 0:
                    info_parts.append(f"PRIORITY={v.priority_score:.2f}")
                info = ";".join(info_parts) if info_parts else "."

                fh.write(
                    f"{v.chromosome}\t{v.position}\t.\t"
                    f"{v.ref}\t{v.alt}\t.\t.\t{info}\n"
                )

        logger.info(f"Wrote {len(collection)} variants to {filepath}")

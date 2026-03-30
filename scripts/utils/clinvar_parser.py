"""
ClinVar VCF parser with star-level stratification.

Parses ClinVar VCF (GRCh38) to extract variants for validation,
filtering by gene, significance, and review status.

Star-level mapping (CLNREVSTAT):
  0 stars: no_assertion, no_criteria
  1 star:  criteria_provided,_single_submitter / criteria_provided,_conflicting
  2 stars: criteria_provided,_multiple_submitters,_no_conflicts
  3 stars: reviewed_by_expert_panel
  4 stars: practice_guideline
"""

import gzip
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from .sequence_utils import Variant

# Default ClinVar VCF path on HPC (shared data)
DEFAULT_CLINVAR_VCF = "/path/to/your/scratch/evo2/data/clinvar/clinvar.vcf.gz"
# Note: Can also import from utils.config: CLINVAR_VCF


# =============================================================================
# Star-level mapping
# =============================================================================

STAR_MAPPING = {
    "no_assertion_criteria_provided": 0,
    "no_assertion_provided": 0,
    "no_classification_provided": 0,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "reviewed_by_expert_panel": 3,
    "practice_guideline": 4,
}

# Clinical significance categories
PATHOGENIC_TERMS = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
BENIGN_TERMS = {"Benign", "Likely_benign", "Benign/Likely_benign"}
VUS_TERMS = {"Uncertain_significance"}

SIGNIFICANCE_MAP = {}
for term in PATHOGENIC_TERMS:
    SIGNIFICANCE_MAP[term] = "P/LP"
for term in BENIGN_TERMS:
    SIGNIFICANCE_MAP[term] = "B/LB"
for term in VUS_TERMS:
    SIGNIFICANCE_MAP[term] = "VUS"


# =============================================================================
# Parser
# =============================================================================

@dataclass
class ClinVarEntry:
    """A parsed ClinVar variant."""
    chrom: str
    pos: int           # 0-based
    ref: str
    alt: str
    clinvar_id: str
    significance: str  # "P/LP", "B/LB", "VUS", or raw
    stars: int
    review_status: str
    gene_symbol: str
    molecular_consequence: str = ""
    origin: str = ""   # germline, somatic, etc.
    raw_significance: str = ""
    submission_date: str = ""  # for temporal splits


def _parse_info(info_str: str) -> Dict[str, str]:
    """Parse VCF INFO field into dict."""
    result = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            result[key] = val
        else:
            result[item] = "True"
    return result


def parse_clinvar_vcf(
    vcf_path: str = DEFAULT_CLINVAR_VCF,
    genes: Optional[Set[str]] = None,
    min_stars: int = 0,
    significance_filter: Optional[Set[str]] = None,
    chrom_filter: Optional[Set[str]] = None,
) -> List[ClinVarEntry]:
    """
    Parse ClinVar VCF and extract variants.

    Args:
        vcf_path: Path to ClinVar VCF (can be .vcf or .vcf.gz)
        genes: Set of gene symbols to filter (None = all)
        min_stars: Minimum review star level (0-4)
        significance_filter: Set of significance categories to include
                           e.g. {"P/LP", "B/LB"} for validation, {"VUS"} for prioritization
        chrom_filter: Set of chromosomes to include (e.g. {"chr17", "chr22"})

    Returns:
        List of ClinVarEntry objects
    """
    entries = []
    path = Path(vcf_path)

    opener = gzip.open if path.suffix == ".gz" else open
    mode = "rt" if path.suffix == ".gz" else "r"

    with opener(vcf_path, mode) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 8:
                continue

            chrom_raw = fields[0]
            # ClinVar uses "1", "2", etc.; normalize to "chr1", "chr2"
            chrom = chrom_raw if chrom_raw.startswith("chr") else f"chr{chrom_raw}"

            if chrom_filter and chrom not in chrom_filter:
                continue

            pos = int(fields[1]) - 1  # VCF is 1-based → convert to 0-based
            ref = fields[3]
            alt = fields[4]
            info = _parse_info(fields[7])

            # Gene filter
            gene_symbol = info.get("GENEINFO", "").split(":")[0] if "GENEINFO" in info else ""
            if genes and gene_symbol not in genes:
                continue

            # Significance
            raw_sig = info.get("CLNSIG", "")
            sig = SIGNIFICANCE_MAP.get(raw_sig, raw_sig)

            if significance_filter and sig not in significance_filter:
                continue

            # Star level
            review_status = info.get("CLNREVSTAT", "")
            # CLNREVSTAT can have multiple entries separated by commas that are part of the value
            stars = STAR_MAPPING.get(review_status, 0)

            if stars < min_stars:
                continue

            # Molecular consequence
            mc = info.get("MC", "")

            # Origin
            origin = info.get("ORIGIN", "")

            entry = ClinVarEntry(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                clinvar_id=fields[2] if len(fields) > 2 else "",
                significance=sig,
                stars=stars,
                review_status=review_status,
                gene_symbol=gene_symbol,
                molecular_consequence=mc,
                origin=origin,
                raw_significance=raw_sig,
            )
            entries.append(entry)

    return entries


def clinvar_to_variants(entries: List[ClinVarEntry]) -> List[Variant]:
    """Convert ClinVarEntry objects to Variant objects for scoring."""
    variants = []
    for e in entries:
        v = Variant(
            chrom=e.chrom,
            pos=e.pos,
            ref=e.ref,
            alt=e.alt,
            variant_id=e.clinvar_id,
            gene=e.gene_symbol,
            clinvar_class=e.significance,
            clinvar_stars=e.stars,
        )
        variants.append(v)
    return variants


def get_validation_sets(
    vcf_path: str = DEFAULT_CLINVAR_VCF,
    genes: Optional[Set[str]] = None,
) -> Dict[str, Dict[str, List[Variant]]]:
    """
    Get pre-built validation sets per gene.

    Returns dict: {gene_symbol: {"P/LP": [...], "B/LB": [...], "VUS": [...]}}
    Each sub-dict has variants stratified by significance.
    Only ≥2-star variants are in P/LP and B/LB (primary validation).
    """
    # Parse all relevant entries
    entries = parse_clinvar_vcf(vcf_path, genes=genes, min_stars=0)

    result: Dict[str, Dict[str, List[Variant]]] = defaultdict(
        lambda: defaultdict(list)
    )

    for e in entries:
        v = Variant(
            chrom=e.chrom, pos=e.pos, ref=e.ref, alt=e.alt,
            variant_id=e.clinvar_id, gene=e.gene_symbol,
            clinvar_class=e.significance, clinvar_stars=e.stars,
        )

        if e.significance in ("P/LP", "B/LB"):
            if e.stars >= 2:
                result[e.gene_symbol][e.significance].append(v)
            # Store all stars (including 1-star) for sensitivity analysis
            result[e.gene_symbol][f"{e.significance}_all_stars"].append(v)
        elif e.significance == "VUS":
            result[e.gene_symbol]["VUS"].append(v)

    return dict(result)


def summarize_clinvar(
    vcf_path: str = DEFAULT_CLINVAR_VCF,
    genes: Optional[Set[str]] = None,
) -> str:
    """Print summary of ClinVar variants per gene and star level."""
    entries = parse_clinvar_vcf(vcf_path, genes=genes)

    # Count by gene, significance, stars
    counts: Dict[str, Dict[str, Dict[int, int]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(int))
    )

    for e in entries:
        counts[e.gene_symbol][e.significance][e.stars] += 1

    lines = [f"ClinVar Summary (source: {vcf_path})"]
    lines.append("=" * 70)

    for gene in sorted(counts.keys()):
        lines.append(f"\n{gene}:")
        for sig in ["P/LP", "B/LB", "VUS"]:
            if sig in counts[gene]:
                star_counts = counts[gene][sig]
                total = sum(star_counts.values())
                stars_str = ", ".join(
                    f"{s}★:{c}" for s, c in sorted(star_counts.items())
                )
                lines.append(f"  {sig:>5}: {total:>5} ({stars_str})")

    return "\n".join(lines)

#!/usr/bin/env python3
"""
Non-coding variant scoring and validation (Phase 5).

Performs:
  A) TERT promoter MPRA validation (Kircher 2019)
  B) ENCODE cCRE enrichment permutation test
  C) Non-coding tool comparison (CADD)

Produces: Figure 6 data

Usage:
  python 05_score_noncoding.py --gene TERT --mode mpra    # TERT MPRA validation
  python 05_score_noncoding.py --gene all  --mode encode   # ENCODE enrichment
  python 05_score_noncoding.py --gene TERT --mode all      # Both
"""

import argparse
import csv
import json
import logging
import random
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from utils.config import (
    GENOME_PATH,
    MPRA_DIR,
    RESULTS_DIR,
    SHARED_DATA_DIR,
)
from utils.gene_coordinates import GENES, TERT, get_gene
from utils.sequence_utils import (
    GenomeAccessor,
    ScoringCheckpoint,
    Variant,
)
from utils.hgvs_mapping import (
    map_tert_mpra,
    parse_hgvs_noncoding,
)

try:
    from scipy import stats as scipy_stats
except ImportError:
    scipy_stats = None

try:
    from sklearn.metrics import roc_auc_score
except ImportError:
    roc_auc_score = None

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)


# =============================================================================
# TERT MPRA Validation
# =============================================================================

TERT_MPRA_GBM_FILE = MPRA_DIR / "tert" / "tert_gbm_scores.csv"
TERT_MPRA_HEK_FILE = MPRA_DIR / "tert" / "tert_hek_scores.csv"

# TERT MPRA reference region (Kircher 2019)
# MaveDB URN: urn:mavedb:00000031-b-1 (GBM primary), urn:mavedb:00000031-a-1 (HEK)
# The MPRA tests a 259bp region of the TERT core promoter.
# Target sequence is on the FORWARD strand (verified by sequence matching).
# Position n.1 = 0-based genomic position 1,294,988 (forward strand, plus direction).
# Position n.259 = 0-based genomic position 1,295,246.
# IMPORTANT: Although TERT is on the minus strand, the MPRA construct uses
# the forward strand. Pass strand="+" to noncoding_hgvs_to_genomic.
TERT_MPRA_REFERENCE_START = 1294988  # 0-based, forward strand n.1


def load_mpra_data() -> List[dict]:
    """Load TERT MPRA scores from CSV."""
    rows = []
    with open(TERT_MPRA_GBM_FILE) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("score") and row["score"] != "NA":
                rows.append(row)
    return rows


def validate_tert_mpra(
    window_size: int = 8192,
) -> Optional[dict]:
    """
    Validate Evo2 scores against TERT MPRA functional data.

    Loads Evo2 scoring results for TERT promoter variants,
    then computes Spearman correlation with MPRA expression scores.
    """
    if not TERT_MPRA_GBM_FILE.exists():
        log.warning("TERT MPRA file not found: {}".format(TERT_MPRA_GBM_FILE))
        return None

    # Load MPRA data and map to genomic coordinates
    mpra_rows = load_mpra_data()
    log.info("MPRA rows: {}".format(len(mpra_rows)))

    mpra_variants = map_tert_mpra(mpra_rows, TERT_MPRA_REFERENCE_START, TERT)
    log.info("MPRA variants mapped: {}".format(len(mpra_variants)))

    if not mpra_variants:
        log.error("No MPRA variants mapped")
        return None

    # Build MPRA score lookup: variant_key -> mpra_score
    mpra_lookup = {}
    for gv in mpra_variants:
        key = "{}:{}:{}>{}".format(gv.chrom, gv.pos, gv.ref, gv.alt)
        mpra_lookup[key] = gv.dms_score

    # Load Evo2 results for TERT
    ckpt_dir = RESULTS_DIR / "TERT" / "w{}".format(window_size)
    checkpoint = ScoringCheckpoint(str(ckpt_dir))
    evo2_results = checkpoint.load_results()

    if not evo2_results:
        log.warning("No Evo2 results for TERT at w{}".format(window_size))
        return None

    # Match Evo2 results with MPRA scores
    matched_evo2 = []
    matched_mpra = []
    matched_keys = []

    for r in evo2_results:
        key = "{}:{}:{}>{}".format(r["chrom"], r["pos"], r["ref"], r["alt"])
        if key in mpra_lookup:
            matched_evo2.append(-r["delta"])  # Negate: higher = more damaging
            matched_mpra.append(mpra_lookup[key])
            matched_keys.append(key)

    log.info("Matched MPRA-Evo2 variants: {}".format(len(matched_evo2)))

    if len(matched_evo2) < 10:
        log.warning("Too few matched variants for correlation")
        return None

    evo2_arr = np.array(matched_evo2)
    mpra_arr = np.array(matched_mpra)

    result = {
        "n_mpra_total": len(mpra_rows),
        "n_mpra_mapped": len(mpra_variants),
        "n_matched": len(matched_evo2),
        "evo2_mean": float(np.mean(evo2_arr)),
        "evo2_std": float(np.std(evo2_arr)),
        "mpra_mean": float(np.mean(mpra_arr)),
        "mpra_std": float(np.std(mpra_arr)),
    }

    if scipy_stats:
        rho, p = scipy_stats.spearmanr(evo2_arr, mpra_arr)
        result["spearman_rho"] = float(rho)
        result["spearman_p"] = float(p)
        log.info("Spearman rho = {:.4f} (p = {:.2e})".format(rho, p))

        pearson_r, pearson_p = scipy_stats.pearsonr(evo2_arr, mpra_arr)
        result["pearson_r"] = float(pearson_r)
        result["pearson_p"] = float(pearson_p)
        log.info("Pearson r = {:.4f} (p = {:.2e})".format(pearson_r, pearson_p))

    # CADD comparison on same MPRA variants
    cadd_results = _compare_cadd_mpra(matched_keys, mpra_arr)
    if cadd_results:
        result["cadd_comparison"] = cadd_results

    # Save per-variant data for scatter plot (Figure 6A)
    variant_data = []
    for key, e, m in zip(matched_keys, matched_evo2, matched_mpra):
        variant_data.append({
            "variant_key": key,
            "evo2_score": float(e),
            "mpra_score": float(m),
        })
    result["variant_data"] = variant_data

    return result


def _compare_cadd_mpra(
    variant_keys: List[str],
    mpra_scores: np.ndarray,
) -> Optional[dict]:
    """Compare CADD scores against MPRA on the same variants."""
    try:
        import pysam
    except ImportError:
        return None

    cadd_path = str(SHARED_DATA_DIR / "benchmarks" / "cadd" / "whole_genome_SNVs.tsv.gz")
    if not Path(cadd_path).exists():
        return None

    # Parse variant keys to get positions
    variant_set = set(variant_keys)
    cadd_scores = {}

    try:
        tbx = pysam.TabixFile(cadd_path)
        # TERT promoter region on chr5
        for row in tbx.fetch("5", 1295100, 1295400):
            fields = row.split("\t")
            chrom = "chr" + fields[0]
            pos = int(fields[1]) - 1
            key = "{}:{}:{}>{}".format(chrom, pos, fields[2], fields[3])
            if key in variant_set:
                cadd_scores[key] = float(fields[5])  # PHRED score
        tbx.close()
    except Exception as e:
        log.warning("CADD MPRA comparison failed: {}".format(e))
        return None

    if len(cadd_scores) < 10:
        return None

    # Match CADD with MPRA
    cadd_matched = []
    mpra_matched = []
    for i, key in enumerate(variant_keys):
        if key in cadd_scores:
            cadd_matched.append(cadd_scores[key])
            mpra_matched.append(mpra_scores[i])

    if len(cadd_matched) < 10 or not scipy_stats:
        return None

    rho, p = scipy_stats.spearmanr(cadd_matched, mpra_matched)
    return {
        "n_matched": len(cadd_matched),
        "spearman_rho": float(rho),
        "spearman_p": float(p),
    }


# =============================================================================
# ENCODE Enrichment Test
# =============================================================================

ENCODE_DIR = SHARED_DATA_DIR / "encode"
CCRE_CLASSES = ["PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"]


def load_ccre_intervals(ccre_class: str) -> List[Tuple[str, int, int]]:
    """Load ENCODE cCRE intervals for a class."""
    bed_file = ENCODE_DIR / "cCRE_{}.bed".format(ccre_class)
    if not bed_file.exists():
        return []
    intervals = []
    with open(bed_file) as f:
        for line in f:
            fields = line.strip().split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            intervals.append((chrom, start, end))
    return intervals


def positions_in_intervals(
    positions: List[Tuple[str, int]],
    intervals: List[Tuple[str, int, int]],
) -> int:
    """Count how many positions fall within any interval."""
    # Build interval lookup by chrom
    by_chrom = defaultdict(list)
    for chrom, start, end in intervals:
        by_chrom[chrom].append((start, end))

    # Sort intervals for binary search
    for chrom in by_chrom:
        by_chrom[chrom].sort()

    count = 0
    for chrom, pos in positions:
        if chrom not in by_chrom:
            continue
        chrom_intervals = by_chrom[chrom]
        # Binary search for overlapping interval
        lo, hi = 0, len(chrom_intervals) - 1
        while lo <= hi:
            mid = (lo + hi) // 2
            start, end = chrom_intervals[mid]
            if pos < start:
                hi = mid - 1
            elif pos >= end:
                lo = mid + 1
            else:
                count += 1
                break

    return count


def compute_encode_enrichment(
    gene_symbol: str,
    window_size: int = 8192,
    top_fraction: float = 0.10,
    n_permutations: int = 1000,
) -> Optional[dict]:
    """
    Test enrichment of Evo2-constrained positions in ENCODE cCREs.

    Method:
    1. Load all scored positions for a gene
    2. Identify top X% most constrained (highest -delta, i.e. most negative delta)
    3. Count overlap with each ENCODE cCRE class
    4. Permutation test: randomly sample same number of positions, compute overlap
    5. Report fold-enrichment and p-value
    """
    gene = get_gene(gene_symbol)

    # Load Evo2 results
    ckpt_dir = RESULTS_DIR / gene_symbol / "w{}".format(window_size)
    checkpoint = ScoringCheckpoint(str(ckpt_dir))
    results = checkpoint.load_results()

    if not results:
        log.warning("No Evo2 results for {}".format(gene_symbol))
        return None

    # Extract positions and scores
    all_positions = []
    all_scores = []
    for r in results:
        chrom = r["chrom"]
        pos = r["pos"]
        delta = r["delta"]
        all_positions.append((chrom, pos))
        all_scores.append(-delta)  # Higher = more constrained

    all_scores = np.array(all_scores)
    n_total = len(all_positions)

    # Top constrained positions
    n_top = max(1, int(n_total * top_fraction))
    top_idx = np.argsort(all_scores)[-n_top:]
    top_positions = [all_positions[i] for i in top_idx]

    log.info("  {}: {} total positions, {} top {}%".format(
        gene_symbol, n_total, n_top, int(top_fraction * 100)))

    enrichment_results = {}

    for ccre_class in CCRE_CLASSES:
        intervals = load_ccre_intervals(ccre_class)
        if not intervals:
            continue

        # Observed overlap
        obs_count = positions_in_intervals(top_positions, intervals)

        # Permutation test
        perm_counts = []
        for _ in range(n_permutations):
            perm_idx = random.sample(range(n_total), n_top)
            perm_positions = [all_positions[i] for i in perm_idx]
            perm_count = positions_in_intervals(perm_positions, intervals)
            perm_counts.append(perm_count)

        perm_counts = np.array(perm_counts)
        mean_perm = np.mean(perm_counts)
        fold_enrichment = obs_count / mean_perm if mean_perm > 0 else float("inf")
        p_value = np.mean(perm_counts >= obs_count)

        enrichment_results[ccre_class] = {
            "observed": int(obs_count),
            "expected_mean": float(mean_perm),
            "expected_std": float(np.std(perm_counts)),
            "fold_enrichment": float(fold_enrichment),
            "p_value": float(p_value),
            "n_permutations": n_permutations,
        }

        log.info("    {}: obs={}, exp={:.1f}, fold={:.2f}, p={:.4f}".format(
            ccre_class, obs_count, mean_perm, fold_enrichment, p_value))

    if not enrichment_results:
        return None

    return {
        "gene": gene_symbol,
        "n_positions": n_total,
        "n_top": n_top,
        "top_fraction": top_fraction,
        "enrichment": enrichment_results,
    }


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Non-coding scoring and validation")
    parser.add_argument("--gene", type=str, default="TERT")
    parser.add_argument("--window-size", type=int, default=8192)
    parser.add_argument("--mode", type=str, default="all",
                        choices=["mpra", "encode", "all"])
    parser.add_argument("--top-fraction", type=float, default=0.10,
                        help="Top fraction for ENCODE enrichment")
    parser.add_argument("--n-permutations", type=int, default=1000)
    args = parser.parse_args()

    out_dir = RESULTS_DIR / "noncoding"
    out_dir.mkdir(parents=True, exist_ok=True)

    all_results = {}

    # A) TERT MPRA validation
    if args.mode in ("mpra", "all"):
        log.info("=" * 60)
        log.info("TERT MPRA Validation")
        log.info("=" * 60)

        mpra_result = validate_tert_mpra(args.window_size)
        if mpra_result:
            all_results["tert_mpra"] = mpra_result

            # Save per-variant scatter data
            if "variant_data" in mpra_result:
                import pandas as pd
                df = pd.DataFrame(mpra_result["variant_data"])
                df.to_csv(out_dir / "tert_mpra_scatter.csv", index=False)

            # Save summary (without bulky variant_data)
            summary = {k: v for k, v in mpra_result.items() if k != "variant_data"}
            with open(out_dir / "tert_mpra_summary.json", "w") as f:
                json.dump(summary, f, indent=2, default=str)

            log.info("MPRA results saved to {}".format(out_dir))
        else:
            log.warning("MPRA validation failed")

    # B) ENCODE enrichment
    if args.mode in ("encode", "all"):
        log.info("\n" + "=" * 60)
        log.info("ENCODE cCRE Enrichment")
        log.info("=" * 60)

        all_genes = [g.symbol for g in GENES.values()]
        genes = all_genes if args.gene.lower() == "all" else [args.gene.upper()]

        encode_results = []
        for gene_sym in genes:
            result = compute_encode_enrichment(
                gene_sym,
                window_size=args.window_size,
                top_fraction=args.top_fraction,
                n_permutations=args.n_permutations,
            )
            if result:
                encode_results.append(result)

        if encode_results:
            all_results["encode_enrichment"] = encode_results

            with open(out_dir / "encode_enrichment.json", "w") as f:
                json.dump(encode_results, f, indent=2, default=str)

            # Summary CSV for plotting (Figure 6B)
            rows = []
            for er in encode_results:
                for ccre_class, enrich in er["enrichment"].items():
                    rows.append({
                        "gene": er["gene"],
                        "ccre_class": ccre_class,
                        "fold_enrichment": enrich["fold_enrichment"],
                        "p_value": enrich["p_value"],
                        "observed": enrich["observed"],
                        "expected_mean": enrich["expected_mean"],
                    })
            import pandas as pd
            pd.DataFrame(rows).to_csv(out_dir / "encode_enrichment.csv", index=False)

            log.info("ENCODE results saved to {}".format(out_dir))

    # Save combined results
    with open(out_dir / "noncoding_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    log.info("\nAll results saved to {}".format(out_dir))


if __name__ == "__main__":
    main()

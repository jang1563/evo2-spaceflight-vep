#!/usr/bin/env python3
"""
Window size ablation study (Phase 1).

Score BRCA1, TP53, CHEK2 ClinVar P/LP and B/LB variants at 5 window sizes:
  4K, 8K, 16K, 32K, 64K

Models:
  - 4K, 8K: evo2_7b (8K native context)
  - 16K-64K: evo2_7b_262k (262K context)

For each window size, compute:
  - AUROC against ClinVar P/LP vs B/LB
  - Score stability (Pearson r between adjacent sizes)
  - Compute time per variant

Produces: Figure 2 data

Usage:
  python 00_window_ablation.py --gene BRCA1 --window-size 8192
  python 00_window_ablation.py --gene all --window-size all
"""

import argparse
import json
import logging
import time
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from utils.config import CLINVAR_VCF, GENOME_PATH, RESULTS_DIR as _RESULTS_DIR
from utils.gene_coordinates import get_gene
from utils.sequence_utils import (
    GenomeAccessor,
    ScoringCheckpoint,
    ScoringResult,
    Variant,
    build_scoring_window,
)
from utils.clinvar_parser import parse_clinvar_vcf, clinvar_to_variants
from utils.benchmarking import compute_auroc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

RESULTS_DIR = _RESULTS_DIR / "window_ablation"

WINDOW_SIZES = [4096, 8192, 16384, 32768, 65536]
ABLATION_GENES = ["BRCA1", "TP53", "CHEK2"]

# Model selection based on window size
# evo2_7b: up to 8K native context
# evo2_7b_262k: up to 262K context
MODEL_FOR_WINDOW = {
    4096: "evo2_7b",
    8192: "evo2_7b",
    16384: "evo2_7b_262k",
    32768: "evo2_7b_262k",
    65536: "evo2_7b_262k",
}


def get_clinvar_validation_variants(gene_symbol: str) -> tuple:
    """
    Get ClinVar P/LP and B/LB variants (≥2 star) for validation.

    Returns: (variants, labels) where labels is 1=P/LP, 0=B/LB
    """
    entries_plp = parse_clinvar_vcf(
        str(CLINVAR_VCF),
        genes={gene_symbol},
        min_stars=2,
        significance_filter={"P/LP"},
    )
    entries_blb = parse_clinvar_vcf(
        str(CLINVAR_VCF),
        genes={gene_symbol},
        min_stars=2,
        significance_filter={"B/LB"},
    )

    variants_plp = clinvar_to_variants(entries_plp)
    variants_blb = clinvar_to_variants(entries_blb)

    # Only keep SNVs for clean comparison
    variants_plp = [v for v in variants_plp if v.is_snv]
    variants_blb = [v for v in variants_blb if v.is_snv]

    for v in variants_plp:
        v.gene = gene_symbol
    for v in variants_blb:
        v.gene = gene_symbol

    variants = variants_plp + variants_blb
    labels = np.array([1] * len(variants_plp) + [0] * len(variants_blb))

    return variants, labels


def run_ablation(
    gene_symbol: str,
    window_size: int,
    genome: GenomeAccessor,
):
    """Score ClinVar variants for one gene at one window size."""
    from utils.sequence_utils import ScoringCheckpoint

    log.info(f"\n--- {gene_symbol} @ {window_size}bp ---")

    # Get variants
    variants, labels = get_clinvar_validation_variants(gene_symbol)
    log.info(f"  Variants: {len(variants)} (P/LP: {sum(labels)}, B/LB: {sum(1 - labels)})")

    if len(variants) < 10:
        log.warning(f"  Too few variants for {gene_symbol}, skipping.")
        return None

    # Select model
    model_name = MODEL_FOR_WINDOW[window_size]
    log.info(f"  Model: {model_name}")

    # Checkpoint
    ckpt_dir = RESULTS_DIR / gene_symbol / f"w{window_size}"
    checkpoint = ScoringCheckpoint(str(ckpt_dir))

    # Load model (lazy to avoid loading before needed)
    log.info(f"  Loading model {model_name}...")
    from evo2 import Evo2
    model = Evo2(model_name)

    # Score variants
    times = []

    for i, var in enumerate(variants):
        if checkpoint.is_scored(var):
            continue

        try:
            ref_win, alt_win = build_scoring_window(genome, var, window_size)
        except ValueError as e:
            log.warning(f"  Skip {var.key()}: {e}")
            continue

        t0 = time.time()
        scores = model.score_sequences(
            seqs=[ref_win.sequence, alt_win.sequence],
            batch_size=1,
            prepend_bos=False,
            reduce_method="mean",
            average_reverse_complement=True,
        )
        elapsed = time.time() - t0

        ref_score, alt_score = float(scores[0]), float(scores[1])

        result = ScoringResult(
            variant=var,
            ref_score=ref_score,
            alt_score=alt_score,
            window_size=window_size,
        )
        checkpoint.save_result(result)

        times.append(elapsed)

        if (i + 1) % 50 == 0:
            log.info(f"  Scored {i + 1}/{len(variants)} ({elapsed:.2f}s/var)")

    # Load all results (including cached)
    all_results = checkpoint.load_results()
    scored_keys = {r["variant_key"]: r["delta"] for r in all_results
                   if "variant_key" in r}

    # Match deltas to labels
    final_deltas = []
    final_labels = []
    for var, label in zip(variants, labels):
        key = var.key()
        if key in scored_keys:
            final_deltas.append(scored_keys[key])
            final_labels.append(label)

    if len(final_deltas) < 10:
        log.warning(f"  Too few scored variants ({len(final_deltas)})")
        return None

    final_deltas = np.array(final_deltas)
    final_labels = np.array(final_labels)

    # Compute AUROC (using -delta as pathogenicity score)
    auroc, _, _, _ = compute_auroc(final_labels, -final_deltas)

    avg_time = np.mean(times) if times else 0

    result = {
        "gene": gene_symbol,
        "window_size": window_size,
        "model": model_name,
        "n_variants": len(final_deltas),
        "n_pathogenic": int(sum(final_labels)),
        "n_benign": int(sum(1 - final_labels)),
        "auroc": float(auroc),
        "mean_time_per_variant": float(avg_time),
        "mean_delta_plp": float(np.mean(final_deltas[final_labels == 1])),
        "mean_delta_blb": float(np.mean(final_deltas[final_labels == 0])),
    }

    log.info(f"  AUROC: {auroc:.4f}")
    log.info(f"  Time/variant: {avg_time:.2f}s")

    # Save result
    result_file = RESULTS_DIR / gene_symbol / f"w{window_size}" / "ablation_result.json"
    result_file.parent.mkdir(parents=True, exist_ok=True)
    with open(result_file, 'w') as f:
        json.dump(result, f, indent=2)

    return result


def main():
    parser = argparse.ArgumentParser(description="Window size ablation study")
    parser.add_argument("--gene", type=str, default="all",
                        help="Gene symbol or 'all' (default: all 3 ablation genes)")
    parser.add_argument("--window-size", type=str, default="all",
                        help="Window size or 'all' (default: all 5 sizes)")
    args = parser.parse_args()

    genes = ABLATION_GENES if args.gene == "all" else [args.gene.upper()]
    windows = WINDOW_SIZES if args.window_size == "all" else [int(args.window_size)]

    log.info(f"Window ablation study")
    log.info(f"  Genes: {genes}")
    log.info(f"  Windows: {windows}")

    genome = GenomeAccessor(str(GENOME_PATH))

    all_results = []
    for gene_sym in genes:
        for ws in windows:
            result = run_ablation(gene_sym, ws, genome)
            if result:
                all_results.append(result)

    genome.close()

    # Summary table
    log.info(f"\n{'='*70}")
    log.info("Ablation Summary")
    log.info(f"{'Gene':<8} {'Window':>8} {'Model':<15} {'AUROC':>8} {'Time/var':>10}")
    log.info("-" * 55)
    for r in all_results:
        log.info(
            f"{r['gene']:<8} {r['window_size']:>8} {r['model']:<15} "
            f"{r['auroc']:>8.4f} {r['mean_time_per_variant']:>10.2f}s"
        )

    # Save aggregate results
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_DIR / "ablation_summary.json", 'w') as f:
        json.dump(all_results, f, indent=2)

    log.info(f"\nResults saved to {RESULTS_DIR / 'ablation_summary.json'}")


if __name__ == "__main__":
    main()

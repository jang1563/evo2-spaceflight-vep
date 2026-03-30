#!/usr/bin/env python3
"""
Core Evo2 VEP scoring engine with checkpointing.

Scores variants using Evo2 log-likelihood ratio:
  delta = score(alt_seq) - score(ref_seq)  [mean autoregressive log-prob]
  pathogenicity_priority = -delta  [more negative delta = higher priority]

Features:
- Reference window deduplication (~35% compute savings)
- Checkpointing for resume on interruption
- Reverse complement averaging
- Configurable window size
- Deterministic at bfloat16 (verified)

Usage:
  python 02_score_variants.py --gene BRCA1 --window-size 8192
  python 02_score_variants.py --gene BRCA1 --variants-file /path/to/variants.jsonl
  python 02_score_variants.py --gene all --window-size 8192
"""

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np

# Add parent dir to path for utils
sys.path.insert(0, str(Path(__file__).parent))

from utils.config import GENOME_PATH, RESULTS_DIR
from utils.gene_coordinates import GENES, get_gene
from utils.sequence_utils import (
    GenomeAccessor,
    ScoringCheckpoint,
    ScoringResult,
    Variant,
    build_scoring_window,
    deduplicate_references,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# =============================================================================
# Scoring Engine
# =============================================================================

class Evo2Scorer:
    """Wrapper around Evo2 model for VEP scoring."""

    def __init__(self, model_name: str = "evo2_7b", device: str = "cuda:0"):
        self.model_name = model_name
        self.device = device
        self.model = None
        self._load_model()

    def _load_model(self):
        log.info(f"Loading Evo2 model: {self.model_name}")
        t0 = time.time()

        from evo2 import Evo2
        self.model = Evo2(self.model_name)

        log.info(f"Model loaded in {time.time() - t0:.1f}s")

    def score_sequences(
        self,
        sequences: list,
        batch_size: int = 1,
        average_rc: bool = True,
    ) -> list:
        """
        Score sequences using Evo2.

        Args:
            sequences: List of DNA sequences
            batch_size: Batch size (1 recommended for A40/A100)
            average_rc: Average with reverse complement (True for VEP)

        Returns:
            List of mean log-probability scores
        """
        scores = self.model.score_sequences(
            seqs=sequences,
            batch_size=batch_size,
            prepend_bos=False,
            reduce_method="mean",
            average_reverse_complement=average_rc,
        )
        return scores

    def score_variant(
        self,
        ref_seq: str,
        alt_seq: str,
        average_rc: bool = True,
    ) -> tuple:
        """Score a single variant (ref vs alt)."""
        scores = self.score_sequences(
            [ref_seq, alt_seq],
            average_rc=average_rc,
        )
        return scores[0], scores[1]


# =============================================================================
# Main Scoring Loop
# =============================================================================

def score_variants(
    variants: list,
    genome: GenomeAccessor,
    scorer: Evo2Scorer,
    checkpoint: ScoringCheckpoint,
    window_size: int = 8192,
    average_rc: bool = True,
) -> list:
    """
    Score a list of variants with deduplication and checkpointing.

    Returns list of ScoringResult objects.
    """
    # Filter out already-scored variants
    remaining = [v for v in variants if not checkpoint.is_scored(v)]
    log.info(
        f"Variants: {len(variants)} total, {len(remaining)} remaining "
        f"({len(variants) - len(remaining)} already scored)"
    )

    if not remaining:
        return []

    # Deduplicate reference windows
    ref_groups = deduplicate_references(remaining, window_size)
    log.info(
        f"Reference windows: {len(ref_groups)} unique "
        f"(dedup saves {len(remaining) - len(ref_groups)} ref scores)"
    )

    results = []
    total_time = 0.0
    n_scored = 0

    for group_idx, (ref_key, group_variants) in enumerate(ref_groups.items()):
        # Build reference window (shared across group)
        first_var = group_variants[0]
        try:
            ref_window, _ = build_scoring_window(genome, first_var, window_size)
        except ValueError as e:
            log.warning(f"Skipping group {ref_key}: {e}")
            continue

        # Score reference once
        ref_scores = scorer.score_sequences(
            [ref_window.sequence], average_rc=average_rc
        )
        ref_score = ref_scores[0]

        # Score each alternate in this group
        t0 = time.time()
        for var in group_variants:
            if checkpoint.is_scored(var):
                continue

            try:
                _, alt_window = build_scoring_window(genome, var, window_size)
            except ValueError as e:
                log.warning(f"Skipping {var.key()}: {e}")
                continue

            alt_scores = scorer.score_sequences(
                [alt_window.sequence], average_rc=average_rc
            )
            alt_score = alt_scores[0]

            elapsed = time.time() - t0

            result = ScoringResult(
                variant=var,
                ref_score=ref_score,
                alt_score=alt_score,
                window_size=window_size,
            )
            results.append(result)
            checkpoint.save_result(result)

            n_scored += 1
            total_time += elapsed
            t0 = time.time()

            if n_scored % 100 == 0:
                rate = total_time / n_scored
                eta = rate * (len(remaining) - n_scored)
                log.info(
                    f"  Scored {n_scored}/{len(remaining)} "
                    f"({rate:.2f}s/var, ETA {eta/3600:.1f}h) "
                    f"delta={result.delta:.6f}"
                )

    log.info(
        f"Scoring complete: {n_scored} variants in {total_time:.0f}s "
        f"({total_time/max(n_scored,1):.2f}s/var)"
    )
    return results


# =============================================================================
# Variant Loading
# =============================================================================

def load_variants_from_file(filepath: str) -> list:
    """Load variants from JSONL file."""
    variants = []
    with open(filepath) as f:
        for line in f:
            data = json.loads(line)
            variants.append(Variant(
                chrom=data["chrom"],
                pos=data["pos"],
                ref=data["ref"],
                alt=data["alt"],
                variant_id=data.get("variant_id", ""),
                gene=data.get("gene", ""),
                region_type=data.get("region_type", ""),
                clinvar_class=data.get("clinvar_class", ""),
                clinvar_stars=data.get("clinvar_stars", 0),
                dms_score=data.get("dms_score"),
            ))
    return variants


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Evo2 VEP Scoring Engine")
    parser.add_argument("--gene", type=str, default="all",
                        help="Gene symbol or 'all'")
    parser.add_argument("--window-size", type=int, default=8192,
                        help="Scoring window size in bp")
    parser.add_argument("--variants-file", type=str, default=None,
                        help="Path to variants JSONL file (overrides --gene)")
    parser.add_argument("--model", type=str, default="evo2_7b",
                        help="Evo2 model name")
    parser.add_argument("--no-rc", action="store_true",
                        help="Disable reverse complement averaging")
    parser.add_argument("--checkpoint-dir", type=str, default=None,
                        help="Checkpoint directory (default: results/<gene>/)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show variant counts without scoring")
    args = parser.parse_args()

    log.info(f"Gene: {args.gene}")
    log.info(f"Window size: {args.window_size}")
    log.info(f"Model: {args.model}")
    log.info(f"RC averaging: {not args.no_rc}")

    # Load variants
    if args.variants_file:
        log.info(f"Loading variants from {args.variants_file}")
        variants = load_variants_from_file(args.variants_file)
        genes_to_score = list({v.gene for v in variants if v.gene})
    else:
        genes_to_score = list(GENES.keys()) if args.gene == "all" else [args.gene.upper()]
        variants = []
        for gene_sym in genes_to_score:
            var_file = RESULTS_DIR / gene_sym / "variants.jsonl"
            if var_file.exists():
                gene_vars = load_variants_from_file(str(var_file))
                variants.extend(gene_vars)
                log.info(f"  {gene_sym}: {len(gene_vars)} variants")
            else:
                log.warning(f"  {gene_sym}: no variants file at {var_file}")

    log.info(f"Total variants to score: {len(variants)}")

    if args.dry_run:
        log.info("Dry run — exiting.")
        return

    # Initialize
    genome = GenomeAccessor(str(GENOME_PATH))
    scorer = Evo2Scorer(model_name=args.model)

    # Score per gene with per-gene checkpoints (so downstream scripts find results)
    if args.checkpoint_dir:
        # Single explicit checkpoint dir
        ckpt = ScoringCheckpoint(args.checkpoint_dir)
        log.info(f"Checkpoint dir: {args.checkpoint_dir} ({ckpt.n_scored} previously scored)")
        results = score_variants(
            variants=variants, genome=genome, scorer=scorer,
            checkpoint=ckpt, window_size=args.window_size,
            average_rc=not args.no_rc,
        )
    else:
        # Per-gene checkpoints (default)
        results = []
        for gene_sym in genes_to_score:
            gene_vars = [v for v in variants if v.gene == gene_sym]
            if not gene_vars:
                continue
            ckpt_dir = RESULTS_DIR / gene_sym / f"w{args.window_size}"
            ckpt = ScoringCheckpoint(str(ckpt_dir))
            log.info(f"\n--- {gene_sym}: {len(gene_vars)} variants, "
                     f"checkpoint: {ckpt_dir} ({ckpt.n_scored} prev) ---")
            gene_results = score_variants(
                variants=gene_vars, genome=genome, scorer=scorer,
                checkpoint=ckpt, window_size=args.window_size,
                average_rc=not args.no_rc,
            )
            results.extend(gene_results)

    log.info(f"\nTotal new results: {len(results)}")

    # Summary statistics
    if results:
        deltas = [r.delta for r in results]
        log.info(f"Delta stats: mean={np.mean(deltas):.6f}, "
                 f"std={np.std(deltas):.6f}, "
                 f"min={np.min(deltas):.6f}, "
                 f"max={np.max(deltas):.6f}")

    genome.close()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Phase 0: 100-variant throughput benchmark.

Measures actual Evo2 scoring speed on GPU:
- Model load time
- Per-variant scoring time (with and without RC averaging)
- Memory usage
- Verifies deterministic scoring at bfloat16

Usage: python benchmark_throughput.py
"""

import json
import logging
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from utils.config import GENOME_PATH, RESULTS_DIR
from utils.gene_coordinates import get_gene
from utils.sequence_utils import (
    GenomeAccessor,
    Variant,
    build_scoring_window,
    generate_all_snvs,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)


def main():
    import torch

    log.info("=" * 60)
    log.info("Evo2 Throughput Benchmark")
    log.info("=" * 60)

    # GPU info
    if torch.cuda.is_available():
        gpu_name = torch.cuda.get_device_name(0)
        gpu_mem = torch.cuda.get_device_properties(0).total_memory / 1e9
        log.info(f"GPU: {gpu_name} ({gpu_mem:.1f} GB)")
    else:
        log.error("No GPU available!")
        return

    # Load genome
    genome = GenomeAccessor(str(GENOME_PATH))

    # Generate 100 test variants from TP53 coding region
    tp53 = get_gene("TP53")
    first_coding = tp53.coding_exons[0]
    ex_start = max(first_coding.start, tp53.cds_start)
    ex_end = min(first_coding.end, tp53.cds_end)

    all_snvs = generate_all_snvs(genome, tp53.chrom, ex_start, ex_start + 40, "TP53")
    test_variants = all_snvs[:100]
    log.info(f"Test variants: {len(test_variants)} SNVs from TP53")

    # Build windows
    log.info("Building scoring windows...")
    windows_8k = []
    for var in test_variants:
        ref_win, alt_win = build_scoring_window(genome, var, window_size=8192)
        windows_8k.append((ref_win, alt_win))
    log.info(f"Built {len(windows_8k)} window pairs")

    # Load model
    log.info("\nLoading evo2_7b...")
    t0 = time.time()
    from evo2 import Evo2
    model = Evo2("evo2_7b")
    load_time = time.time() - t0
    log.info(f"Model load time: {load_time:.1f}s")
    log.info(f"GPU memory after load: {torch.cuda.memory_allocated()/1e9:.2f} GB")

    # =========================================================================
    # Benchmark 1: Scoring speed WITHOUT RC averaging
    # =========================================================================
    log.info("\n--- Benchmark 1: No RC averaging ---")
    times_no_rc = []
    deltas_no_rc = []

    for i, (ref_win, alt_win) in enumerate(windows_8k):
        t0 = time.time()
        scores = model.score_sequences(
            seqs=[ref_win.sequence, alt_win.sequence],
            batch_size=1,
            prepend_bos=False,
            reduce_method="mean",
            average_reverse_complement=False,
        )
        elapsed = time.time() - t0
        times_no_rc.append(elapsed)
        deltas_no_rc.append(scores[1] - scores[0])

        if (i + 1) % 20 == 0:
            log.info(f"  {i + 1}/100 ({elapsed:.3f}s)")

    log.info(f"  Mean time/variant: {np.mean(times_no_rc):.3f}s")
    log.info(f"  Std time/variant: {np.std(times_no_rc):.3f}s")
    log.info(f"  Total time: {sum(times_no_rc):.1f}s")

    # =========================================================================
    # Benchmark 2: Scoring speed WITH RC averaging
    # =========================================================================
    log.info("\n--- Benchmark 2: With RC averaging ---")
    times_rc = []
    deltas_rc = []

    for i, (ref_win, alt_win) in enumerate(windows_8k):
        t0 = time.time()
        scores = model.score_sequences(
            seqs=[ref_win.sequence, alt_win.sequence],
            batch_size=1,
            prepend_bos=False,
            reduce_method="mean",
            average_reverse_complement=True,
        )
        elapsed = time.time() - t0
        times_rc.append(elapsed)
        deltas_rc.append(scores[1] - scores[0])

        if (i + 1) % 20 == 0:
            log.info(f"  {i + 1}/100 ({elapsed:.3f}s)")

    log.info(f"  Mean time/variant: {np.mean(times_rc):.3f}s")
    log.info(f"  Std time/variant: {np.std(times_rc):.3f}s")
    log.info(f"  Total time: {sum(times_rc):.1f}s")
    log.info(f"  RC overhead: {np.mean(times_rc)/np.mean(times_no_rc):.2f}x")

    # =========================================================================
    # Benchmark 3: Deterministic check at bfloat16
    # =========================================================================
    log.info("\n--- Benchmark 3: Deterministic check ---")
    # Score first 10 variants twice, verify identical
    n_check = 10
    all_identical = True
    for i in range(n_check):
        ref_win, alt_win = windows_8k[i]
        scores1 = model.score_sequences(
            seqs=[ref_win.sequence, alt_win.sequence],
            batch_size=1, prepend_bos=False,
            reduce_method="mean", average_reverse_complement=True,
        )
        scores2 = model.score_sequences(
            seqs=[ref_win.sequence, alt_win.sequence],
            batch_size=1, prepend_bos=False,
            reduce_method="mean", average_reverse_complement=True,
        )
        if scores1[0] != scores2[0] or scores1[1] != scores2[1]:
            log.warning(f"  Variant {i}: NOT deterministic!")
            log.warning(f"    Run1: {scores1[0]:.10f}, {scores1[1]:.10f}")
            log.warning(f"    Run2: {scores2[0]:.10f}, {scores2[1]:.10f}")
            all_identical = False
        else:
            log.info(f"  Variant {i}: deterministic (ref={scores1[0]:.10f})")

    if all_identical:
        log.info("  DETERMINISTIC: All 10 variants produced identical scores")
    else:
        log.warning("  NON-DETERMINISTIC: Some variants produced different scores!")

    # =========================================================================
    # Summary
    # =========================================================================
    log.info("\n" + "=" * 60)
    log.info("BENCHMARK SUMMARY")
    log.info("=" * 60)
    log.info(f"GPU: {gpu_name}")
    log.info(f"Model: evo2_7b")
    log.info(f"Window size: 8192 bp")
    log.info(f"Variants tested: {len(test_variants)}")
    log.info(f"Model load time: {load_time:.1f}s")
    log.info(f"Time/variant (no RC): {np.mean(times_no_rc):.3f}s ± {np.std(times_no_rc):.3f}s")
    log.info(f"Time/variant (with RC): {np.mean(times_rc):.3f}s ± {np.std(times_rc):.3f}s")
    log.info(f"Deterministic: {all_identical}")
    log.info(f"GPU memory: {torch.cuda.memory_allocated()/1e9:.2f} GB allocated")
    log.info(f"Delta range: [{min(deltas_rc):.6f}, {max(deltas_rc):.6f}]")
    log.info(f"Delta mean: {np.mean(deltas_rc):.6f}")

    # Extrapolate total project compute
    time_per_var_rc = np.mean(times_rc)
    # Rough estimate: ~200K total variants across all genes
    est_total_vars = 200000
    est_hours = est_total_vars * time_per_var_rc / 3600
    log.info(f"\nEstimated total project time @ {time_per_var_rc:.3f}s/var:")
    log.info(f"  ~{est_total_vars:,} variants → ~{est_hours:.0f} GPU-hours")

    # Save results
    results = {
        "gpu": gpu_name,
        "gpu_memory_gb": gpu_mem,
        "model": "evo2_7b",
        "window_size": 8192,
        "n_variants": len(test_variants),
        "model_load_time_s": load_time,
        "mean_time_no_rc_s": float(np.mean(times_no_rc)),
        "std_time_no_rc_s": float(np.std(times_no_rc)),
        "mean_time_rc_s": float(np.mean(times_rc)),
        "std_time_rc_s": float(np.std(times_rc)),
        "rc_overhead": float(np.mean(times_rc) / np.mean(times_no_rc)),
        "deterministic": all_identical,
        "gpu_memory_allocated_gb": float(torch.cuda.memory_allocated() / 1e9),
        "delta_mean": float(np.mean(deltas_rc)),
        "delta_std": float(np.std(deltas_rc)),
        "delta_min": float(min(deltas_rc)),
        "delta_max": float(max(deltas_rc)),
    }

    out_file = RESULTS_DIR / "benchmark_throughput.json"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, 'w') as f:
        json.dump(results, f, indent=2)
    log.info(f"\nResults saved to {out_file}")

    genome.close()


if __name__ == "__main__":
    main()

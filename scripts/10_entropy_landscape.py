#!/usr/bin/env python3
"""
Positional entropy landscape computation (Phase 5 supplement).

Computes per-position entropy across gene regions using Evo2's
positional_entropies function. Overlays with ENCODE cCREs to
discover cryptic constrained non-coding elements.

Produces: entropy landscape data for Figure 5 panels

Usage:
  python 10_entropy_landscape.py --gene TERT --window-size 8192
  python 10_entropy_landscape.py --gene all
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from utils.config import GENOME_PATH, RESULTS_DIR, SHARED_DATA_DIR
from utils.gene_coordinates import GENES, get_gene

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)


# =============================================================================
# Entropy Computation
# =============================================================================

def compute_gene_entropy(
    gene_symbol: str,
    window_size: int = 8192,
    stride: int = 4096,
    padding: int = 5000,
) -> Optional[dict]:
    """
    Compute positional entropy across a gene region using Evo2.

    Tiles the gene region with overlapping windows, computes per-position
    entropy, then merges overlapping positions by averaging.

    Args:
        gene_symbol: Gene symbol
        window_size: Window size for entropy computation
        stride: Step size between windows (< window_size for overlap)
        padding: bp to extend beyond gene boundaries
    """
    gene = get_gene(gene_symbol)
    chrom, region_start, region_end = gene.get_region(padding)
    region_length = region_end - region_start

    log.info("{}: {}:{}-{} ({} bp)".format(
        gene_symbol, chrom, region_start, region_end, region_length))

    # Check for existing results
    out_dir = RESULTS_DIR / gene_symbol / "entropy"
    out_dir.mkdir(parents=True, exist_ok=True)
    result_file = out_dir / "entropy_w{}.npz".format(window_size)

    if result_file.exists():
        log.info("  Loading cached entropy from {}".format(result_file))
        data = np.load(result_file)
        return {
            "gene": gene_symbol,
            "chrom": chrom,
            "start": region_start,
            "end": region_end,
            "positions": data["positions"].tolist(),
            "entropies": data["entropies"].tolist(),
            "window_size": window_size,
        }

    # Import Evo2 and genome accessor (GPU required)
    try:
        from evo2 import Evo2
        from evo2.scoring import positional_entropies
        from utils.sequence_utils import GenomeAccessor
    except ImportError as e:
        log.error("Cannot import Evo2 or pyfaidx: {}".format(e))
        return None

    genome = GenomeAccessor(str(GENOME_PATH))

    log.info("Loading Evo2 model...")
    model_obj = Evo2("evo2_7b")

    # Tile the region with overlapping windows
    n_windows = max(1, (region_length - window_size) // stride + 1)
    log.info("  {} windows (size={}, stride={})".format(
        n_windows, window_size, stride))

    # Accumulate entropy values per position
    entropy_sum = np.zeros(region_length, dtype=np.float64)
    entropy_count = np.zeros(region_length, dtype=np.int32)

    t0 = time.time()
    for i in range(n_windows):
        win_start = region_start + i * stride
        win_end = win_start + window_size

        # Don't go past region end
        if win_end > region_end + padding:
            break

        seq = genome.fetch(chrom, win_start, win_end)
        if len(seq) < window_size:
            continue

        # Compute positional entropies
        entropies = positional_entropies(
            [seq], model_obj.model, model_obj.tokenizer
        )
        ent_arr = np.array(entropies[0])

        # Map to region coordinates
        for j in range(len(ent_arr)):
            pos_in_region = (win_start - region_start) + j
            if 0 <= pos_in_region < region_length:
                entropy_sum[pos_in_region] += ent_arr[j]
                entropy_count[pos_in_region] += 1

        if (i + 1) % 10 == 0:
            elapsed = time.time() - t0
            log.info("  Window {}/{} ({:.1f}s)".format(i + 1, n_windows, elapsed))

    genome.close()

    # Average overlapping positions
    valid = entropy_count > 0
    mean_entropy = np.zeros(region_length)
    mean_entropy[valid] = entropy_sum[valid] / entropy_count[valid]

    positions = np.arange(region_start, region_end)

    # Save as npz for later plotting
    np.savez_compressed(
        result_file,
        positions=positions,
        entropies=mean_entropy[:len(positions)],
        entropy_count=entropy_count[:len(positions)],
    )
    log.info("  Saved to {}".format(result_file))

    elapsed = time.time() - t0
    log.info("  Computed in {:.0f}s ({:.2f}s/window)".format(
        elapsed, elapsed / max(n_windows, 1)))

    return {
        "gene": gene_symbol,
        "chrom": chrom,
        "start": region_start,
        "end": region_end,
        "n_positions": int(np.sum(valid)),
        "mean_entropy": float(np.mean(mean_entropy[valid])),
        "window_size": window_size,
        "n_windows": n_windows,
        "compute_time_s": elapsed,
    }


# =============================================================================
# ENCODE Overlay Analysis
# =============================================================================

def overlay_encode(
    gene_symbol: str,
    window_size: int = 8192,
) -> Optional[dict]:
    """
    Overlay entropy landscape with ENCODE cCREs.

    Computes mean entropy inside vs outside each cCRE class.
    """
    gene = get_gene(gene_symbol)
    chrom, region_start, region_end = gene.get_region(5000)

    # Load entropy
    ent_file = RESULTS_DIR / gene_symbol / "entropy" / "entropy_w{}.npz".format(window_size)
    if not ent_file.exists():
        return None

    data = np.load(ent_file)
    positions = data["positions"]
    entropies = data["entropies"]

    if len(positions) == 0:
        return None

    # Create position -> entropy lookup
    pos_to_ent = dict(zip(positions, entropies))

    encode_dir = SHARED_DATA_DIR / "encode"
    ccre_classes = ["PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"]

    results = {}
    for ccre_class in ccre_classes:
        bed_file = encode_dir / "cCRE_{}.bed".format(ccre_class)
        if not bed_file.exists():
            continue

        # Find cCREs overlapping this gene region
        inside_ent = []
        with open(bed_file) as f:
            for line in f:
                fields = line.strip().split("\t")
                if fields[0] != chrom:
                    continue
                start = int(fields[1])
                end = int(fields[2])
                if end < region_start or start > region_end:
                    continue
                # Collect entropy values inside this cCRE
                for pos in range(max(start, region_start), min(end, region_end)):
                    if pos in pos_to_ent:
                        inside_ent.append(pos_to_ent[pos])

        if not inside_ent:
            continue

        # Outside entropy (all positions not in any cCRE of this class)
        inside_pos = set()
        with open(bed_file) as f:
            for line in f:
                fields = line.strip().split("\t")
                if fields[0] != chrom:
                    continue
                start = int(fields[1])
                end = int(fields[2])
                if end < region_start or start > region_end:
                    continue
                for pos in range(max(start, region_start), min(end, region_end)):
                    inside_pos.add(pos)

        outside_ent = [pos_to_ent[p] for p in positions
                       if p not in inside_pos and p in pos_to_ent]

        if not outside_ent:
            continue

        inside_arr = np.array(inside_ent)
        outside_arr = np.array(outside_ent)

        entry = {
            "n_inside": len(inside_ent),
            "n_outside": len(outside_ent),
            "mean_inside": float(np.mean(inside_arr)),
            "mean_outside": float(np.mean(outside_arr)),
            "ratio": float(np.mean(inside_arr) / np.mean(outside_arr))
            if np.mean(outside_arr) > 0 else float("inf"),
        }

        try:
            from scipy import stats as scipy_stats
            _, p = scipy_stats.mannwhitneyu(inside_arr, outside_arr, alternative="two-sided")
            entry["mannwhitney_p"] = float(p)
        except (ImportError, ValueError):
            pass

        results[ccre_class] = entry

    return results if results else None


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Positional entropy landscape")
    parser.add_argument("--gene", type=str, default="TERT")
    parser.add_argument("--window-size", type=int, default=8192)
    parser.add_argument("--stride", type=int, default=4096)
    parser.add_argument("--padding", type=int, default=5000)
    parser.add_argument("--encode-overlay", action="store_true",
                        help="Run ENCODE overlay analysis on cached entropy data")
    args = parser.parse_args()

    all_symbols = [g.symbol for g in GENES.values()]
    genes = all_symbols if args.gene.lower() == "all" else [args.gene.upper()]

    if args.encode_overlay:
        # Overlay mode: use cached entropy data
        log.info("ENCODE Overlay Analysis")
        for gene_sym in genes:
            result = overlay_encode(gene_sym, args.window_size)
            if result:
                log.info("\n{}: ENCODE cCRE entropy overlay".format(gene_sym))
                for ccre_class, data in result.items():
                    log.info("  {}: inside={:.4f} outside={:.4f} ratio={:.2f} p={:.2e}".format(
                        ccre_class,
                        data["mean_inside"], data["mean_outside"],
                        data["ratio"], data.get("mannwhitney_p", 1),
                    ))

                out_dir = RESULTS_DIR / gene_sym / "entropy"
                with open(out_dir / "encode_overlay.json", "w") as f:
                    json.dump(result, f, indent=2, default=str)
        return

    # Compute entropy (GPU required)
    log.info("Positional Entropy Landscape — genes: {}".format(genes))

    results = []
    for gene_sym in genes:
        result = compute_gene_entropy(
            gene_sym, args.window_size, args.stride, args.padding
        )
        if result:
            results.append(result)

    if results:
        # Summary
        log.info("\n" + "=" * 60)
        log.info("ENTROPY SUMMARY")
        log.info("=" * 60)
        for r in results:
            log.info("  {}: {} positions, mean entropy={:.4f}".format(
                r["gene"], r.get("n_positions", 0), r.get("mean_entropy", 0)))

        # Save summary
        out_dir = RESULTS_DIR / "entropy_landscapes"
        out_dir.mkdir(parents=True, exist_ok=True)
        with open(out_dir / "entropy_summary.json", "w") as f:
            json.dump(results, f, indent=2, default=str)
        log.info("Summary saved to {}".format(out_dir))


if __name__ == "__main__":
    main()

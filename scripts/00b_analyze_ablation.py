#!/usr/bin/env python3
"""
Analyze window size ablation results (Phase 1 post-processing).

Aggregates scoring results from all gene × window combinations,
computes AUROC/AUPRC, score stability, and produces Figure 2 data.

Usage:
  python 00b_analyze_ablation.py              # Analyze all available results
  python 00b_analyze_ablation.py --partial    # Include incomplete runs
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from utils.config import RESULTS_DIR as _RESULTS_DIR
from utils.benchmarking import compute_auroc, compute_auprc

try:
    from scipy import stats as scipy_stats
except ImportError:
    scipy_stats = None

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

RESULTS_DIR = _RESULTS_DIR / "window_ablation"
ABLATION_GENES = ["BRCA1", "TP53", "CHEK2"]
WINDOW_SIZES = [4096, 8192, 16384, 32768, 65536]


def load_checkpoint_results(gene: str, window_size: int) -> list:
    """Load scored variants from checkpoint JSONL file."""
    ckpt_file = RESULTS_DIR / gene / f"w{window_size}" / "results.jsonl"
    if not ckpt_file.exists():
        return []
    results = []
    with open(ckpt_file) as f:
        for line in f:
            results.append(json.loads(line.strip()))
    return results


def compute_metrics(results: list, min_variants: int = 10) -> dict:
    """Compute AUROC and AUPRC from scored variants with ClinVar labels."""
    plp_deltas = [r["delta"] for r in results if r.get("clinvar_class") == "P/LP"]
    blb_deltas = [r["delta"] for r in results if r.get("clinvar_class") == "B/LB"]

    if len(plp_deltas) < 5 or len(blb_deltas) < 5:
        return None

    labels = np.array([1] * len(plp_deltas) + [0] * len(blb_deltas))
    deltas = np.array(plp_deltas + blb_deltas)
    scores = -deltas  # higher = more pathogenic

    auroc, fpr, tpr, _ = compute_auroc(labels, scores)
    auprc, prec, recall, _ = compute_auprc(labels, scores)

    return {
        "n_total": len(results),
        "n_plp": len(plp_deltas),
        "n_blb": len(blb_deltas),
        "auroc": float(auroc),
        "auprc": float(auprc),
        "mean_delta_plp": float(np.mean(plp_deltas)),
        "mean_delta_blb": float(np.mean(blb_deltas)),
        "std_delta_plp": float(np.std(plp_deltas)),
        "std_delta_blb": float(np.std(blb_deltas)),
        "separation": float(np.mean(blb_deltas) - np.mean(plp_deltas)),
    }


def compute_stability(
    results_small: list, results_large: list
) -> dict:
    """
    Compute score stability between adjacent window sizes.

    Matches variants by key and computes Pearson r of deltas.
    """
    if not scipy_stats:
        return {"pearson_r": None, "n_shared": 0}

    # Build delta lookup for each
    deltas_small = {}
    for r in results_small:
        key = r.get("variant_key", f"{r['chrom']}:{r['pos']}:{r['ref']}>{r['alt']}")
        deltas_small[key] = r["delta"]

    deltas_large = {}
    for r in results_large:
        key = r.get("variant_key", f"{r['chrom']}:{r['pos']}:{r['ref']}>{r['alt']}")
        deltas_large[key] = r["delta"]

    # Find shared variants
    shared_keys = set(deltas_small.keys()) & set(deltas_large.keys())
    if len(shared_keys) < 10:
        return {"pearson_r": None, "n_shared": len(shared_keys)}

    x = np.array([deltas_small[k] for k in shared_keys])
    y = np.array([deltas_large[k] for k in shared_keys])

    r, p = scipy_stats.pearsonr(x, y)
    return {
        "pearson_r": float(r),
        "pearson_p": float(p),
        "n_shared": len(shared_keys),
    }


def analyze_star_stratification(results: list) -> dict:
    """Compute AUROC at different ClinVar star levels."""
    star_aurocs = {}
    for min_stars in [1, 2, 3]:
        plp = [r["delta"] for r in results
               if r.get("clinvar_class") == "P/LP"
               and r.get("clinvar_stars", 0) >= min_stars]
        blb = [r["delta"] for r in results
               if r.get("clinvar_class") == "B/LB"
               and r.get("clinvar_stars", 0) >= min_stars]

        if len(plp) >= 5 and len(blb) >= 5:
            labels = np.array([1] * len(plp) + [0] * len(blb))
            scores = -np.array(plp + blb)
            auroc, _, _, _ = compute_auroc(labels, scores)
            star_aurocs[f">=_{min_stars}_star"] = {
                "auroc": float(auroc),
                "n_plp": len(plp),
                "n_blb": len(blb),
            }
    return star_aurocs


def main():
    parser = argparse.ArgumentParser(description="Analyze ablation results")
    parser.add_argument("--partial", action="store_true",
                        help="Include incomplete runs in analysis")
    args = parser.parse_args()

    log.info("Window Size Ablation Analysis")
    log.info(f"Results dir: {RESULTS_DIR}")

    # Collect all results
    all_metrics = []
    all_results_cache = {}  # (gene, ws) -> results

    for gene in ABLATION_GENES:
        for ws in WINDOW_SIZES:
            results = load_checkpoint_results(gene, ws)
            if not results:
                log.info(f"  {gene} @ {ws}bp: no results yet")
                continue

            all_results_cache[(gene, ws)] = results
            metrics = compute_metrics(results)
            if metrics is None:
                log.warning(f"  {gene} @ {ws}bp: too few labeled variants "
                            f"(n={len(results)})")
                continue

            row = {"gene": gene, "window_size": ws, **metrics}
            all_metrics.append(row)

            # Star stratification
            star_strat = analyze_star_stratification(results)
            row["star_stratification"] = star_strat

            log.info(
                f"  {gene} @ {ws:>6}bp: AUROC={metrics['auroc']:.4f}  "
                f"AUPRC={metrics['auprc']:.4f}  "
                f"n={metrics['n_plp']}P/{metrics['n_blb']}B  "
                f"sep={metrics['separation']:.6f}"
            )

    if not all_metrics:
        log.error("No results available for analysis")
        return

    # Convert to DataFrame
    df = pd.DataFrame(all_metrics)

    # Compute stability between adjacent window sizes
    stability_results = []
    for gene in ABLATION_GENES:
        for i in range(len(WINDOW_SIZES) - 1):
            ws_small = WINDOW_SIZES[i]
            ws_large = WINDOW_SIZES[i + 1]

            if (gene, ws_small) not in all_results_cache:
                continue
            if (gene, ws_large) not in all_results_cache:
                continue

            stab = compute_stability(
                all_results_cache[(gene, ws_small)],
                all_results_cache[(gene, ws_large)],
            )
            stab["gene"] = gene
            stab["window_small"] = ws_small
            stab["window_large"] = ws_large
            stability_results.append(stab)

            if stab["pearson_r"] is not None:
                log.info(
                    f"  Stability {gene} {ws_small}->{ws_large}: "
                    f"r={stab['pearson_r']:.4f} (n={stab['n_shared']})"
                )

    # Summary table
    log.info(f"\n{'='*80}")
    log.info("ABLATION SUMMARY")
    log.info(f"{'='*80}")
    log.info(f"{'Gene':<8} {'Window':>8} {'AUROC':>8} {'AUPRC':>8} "
             f"{'P/LP':>6} {'B/LB':>6} {'Separation':>12}")
    log.info("-" * 62)
    for _, row in df.iterrows():
        log.info(
            f"{row['gene']:<8} {row['window_size']:>8} "
            f"{row['auroc']:>8.4f} {row['auprc']:>8.4f} "
            f"{row['n_plp']:>6} {row['n_blb']:>6} "
            f"{row['separation']:>12.6f}"
        )

    # Optimal window selection
    log.info(f"\n{'='*80}")
    log.info("OPTIMAL WINDOW SELECTION")
    log.info(f"{'='*80}")

    for gene in ABLATION_GENES:
        gene_df = df[df["gene"] == gene]
        if gene_df.empty:
            continue
        best_idx = gene_df["auroc"].idxmax()
        best = gene_df.loc[best_idx]
        log.info(
            f"  {gene}: best AUROC={best['auroc']:.4f} "
            f"at window={int(best['window_size'])}bp"
        )

    # Overall recommendation (mean AUROC across genes per window)
    if len(df["gene"].unique()) > 1:
        mean_auroc = df.groupby("window_size")["auroc"].mean()
        best_ws = mean_auroc.idxmax()
        log.info(
            f"\n  OVERALL: best mean AUROC={mean_auroc[best_ws]:.4f} "
            f"at window={best_ws}bp"
        )

    # Save results
    out_dir = RESULTS_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    # Summary JSON
    summary = {
        "metrics": all_metrics,
        "stability": stability_results,
        "recommendation": {
            "per_gene": {},
            "overall": None,
        },
    }

    for gene in ABLATION_GENES:
        gene_df = df[df["gene"] == gene]
        if not gene_df.empty:
            best_idx = gene_df["auroc"].idxmax()
            summary["recommendation"]["per_gene"][gene] = {
                "best_window": int(gene_df.loc[best_idx, "window_size"]),
                "best_auroc": float(gene_df.loc[best_idx, "auroc"]),
            }

    if len(df["gene"].unique()) > 1:
        mean_auroc = df.groupby("window_size")["auroc"].mean()
        best_ws = int(mean_auroc.idxmax())
        summary["recommendation"]["overall"] = {
            "best_window": best_ws,
            "mean_auroc": float(mean_auroc[best_ws]),
        }

    with open(out_dir / "ablation_analysis.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # CSV for easy plotting
    df.to_csv(out_dir / "ablation_metrics.csv", index=False)

    if stability_results:
        pd.DataFrame(stability_results).to_csv(
            out_dir / "ablation_stability.csv", index=False
        )

    log.info(f"\nResults saved to {out_dir}")
    log.info(f"  ablation_analysis.json — full analysis")
    log.info(f"  ablation_metrics.csv — metrics table")
    if stability_results:
        log.info(f"  ablation_stability.csv — stability data")


if __name__ == "__main__":
    main()

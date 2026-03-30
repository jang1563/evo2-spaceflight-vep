#!/usr/bin/env python3
"""
Publication figure generation (Phase 8).

Generates all main and supplementary figures for the manuscript.
Reads pre-computed results from JSON/CSV files produced by upstream scripts.

Main Figures:
  Fig 1 — Study overview (manual/Illustrator)
  Fig 2 — Window size ablation
  Fig 3 — DMS calibration (4 panels)
  Fig 4 — ClinVar validation + tool comparison
  Fig 5 — Per-gene constraint landscapes
  Fig 6 — Non-coding validation (MPRA + ENCODE)
  Fig 7 — Radiation functional impact
  Fig 8 — Spaceflight integration (astronaut variants)

Supplementary:
  S1 — Orthogonality matrix
  S2 — 7B vs 40B comparison
  S3 — Indel analysis
  S4 — Temporal ClinVar split
  S5 — ClinVar star-level sensitivity
  S6 — Gene-specific Pejaver LR calibration curves

Usage:
  python 12_make_figures.py --fig 2      # Single figure
  python 12_make_figures.py --fig all    # All figures
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from utils.config import RESULTS_DIR
from utils.gene_coordinates import GENES, CONTROL_GENES, NOVEL_GENES

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import TwoSlopeNorm
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    import seaborn as sns
    HAS_SNS = True
except ImportError:
    HAS_SNS = False

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

FIG_DIR = RESULTS_DIR / "figures"

# Publication style
COLORS = {
    "Evo2": "#1f77b4",
    "AlphaMissense": "#ff7f0e",
    "CADD": "#2ca02c",
    "REVEL": "#d62728",
    "LINSIGHT": "#9467bd",
    "SpliceAI": "#8c564b",
    "P/LP": "#d62728",
    "B/LB": "#2ca02c",
    "VUS": "#7f7f7f",
    "radiation": "#e41a1c",
    "other": "#999999",
}


def setup_style():
    """Set publication figure style."""
    if not HAS_MPL:
        return
    plt.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.family": "sans-serif",
    })
    if HAS_SNS:
        sns.set_style("whitegrid")


# =============================================================================
# Figure 2: Window Size Ablation
# =============================================================================

def fig2_ablation():
    """Window size ablation: AUROC vs window size + stability heatmap."""
    csv_file = RESULTS_DIR / "window_ablation" / "ablation_metrics.csv"
    if not csv_file.exists():
        log.warning("No ablation data: {}".format(csv_file))
        return

    df = pd.read_csv(csv_file)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: AUROC vs window size
    ax = axes[0]
    for gene in df["gene"].unique():
        gdf = df[df["gene"] == gene].sort_values("window_size")
        ax.plot(gdf["window_size"], gdf["auroc"], "o-", label=gene, markersize=6)

    ax.set_xlabel("Window Size (bp)")
    ax.set_ylabel("AUROC")
    ax.set_xscale("log", base=2)
    ax.set_xticks([4096, 8192, 16384, 32768, 65536])
    ax.set_xticklabels(["4K", "8K", "16K", "32K", "64K"])
    ax.set_ylim(0.5, 1.0)
    ax.legend()
    ax.set_title("A. AUROC vs Window Size")

    # Panel B: Stability heatmap
    stab_file = RESULTS_DIR / "window_ablation" / "ablation_stability.csv"
    ax = axes[1]
    if stab_file.exists():
        stab_df = pd.read_csv(stab_file)
        stab_df = stab_df.dropna(subset=["pearson_r"])
        if not stab_df.empty:
            pivot = stab_df.pivot_table(
                index="gene", columns="window_large",
                values="pearson_r", aggfunc="first"
            )
            if HAS_SNS:
                sns.heatmap(pivot, ax=ax, annot=True, fmt=".3f",
                           cmap="YlOrRd", vmin=0.9, vmax=1.0)
            ax.set_title("B. Score Stability (Pearson r)")
    else:
        ax.text(0.5, 0.5, "Stability data\nnot yet available",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title("B. Score Stability")

    plt.tight_layout()
    out = FIG_DIR / "fig2_ablation.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    log.info("Saved: {}".format(out))


# =============================================================================
# Figure 3: DMS Calibration
# =============================================================================

def fig3_dms_calibration():
    """DMS calibration scatter plots for 4 control genes."""
    genes = ["BRCA1", "TP53", "CHEK2", "DNMT3A"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    for idx, gene in enumerate(genes):
        ax = axes[idx // 2, idx % 2]
        csv_file = RESULTS_DIR / gene / "dms_calibration.csv"

        if csv_file.exists():
            df = pd.read_csv(csv_file)
            ax.scatter(df["dms_score"], df["evo2_delta"],
                      alpha=0.3, s=8, color=COLORS["Evo2"])

            # Add correlation
            if len(df) > 10:
                from scipy import stats
                rho, p = stats.spearmanr(df["dms_score"], df["evo2_delta"])
                ax.text(0.05, 0.95, "rho = {:.3f}".format(rho),
                       transform=ax.transAxes, va="top", fontsize=10)

            ax.set_xlabel("DMS Score")
            ax.set_ylabel("Evo2 Delta")
        else:
            ax.text(0.5, 0.5, "Data not yet\navailable",
                    ha="center", va="center", transform=ax.transAxes)

        ax.set_title("{}) {}".format(chr(65 + idx), gene))

    plt.tight_layout()
    out = FIG_DIR / "fig3_dms_calibration.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    log.info("Saved: {}".format(out))


# =============================================================================
# Figure 4: ClinVar Validation + Tool Comparison
# =============================================================================

def fig4_clinvar_validation():
    """ClinVar validation ROC curves and tool comparison bars."""
    json_file = RESULTS_DIR / "clinvar_validation" / "validation_results.json"
    tool_csv = RESULTS_DIR / "benchmark_comparison" / "tool_aurocs.csv"

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: AUROC bar chart by tool
    ax = axes[0]
    if tool_csv.exists():
        df = pd.read_csv(tool_csv)
        mean_auroc = df.groupby("tool")["auroc"].agg(["mean", "std"]).reset_index()
        mean_auroc = mean_auroc.sort_values("mean", ascending=True)

        colors = [COLORS.get(t, "#7f7f7f") for t in mean_auroc["tool"]]
        ax.barh(mean_auroc["tool"], mean_auroc["mean"],
                xerr=mean_auroc["std"], color=colors, capsize=3)
        ax.set_xlabel("Mean AUROC")
        ax.set_xlim(0.5, 1.0)
    else:
        ax.text(0.5, 0.5, "Data not yet\navailable",
                ha="center", va="center", transform=ax.transAxes)
    ax.set_title("A. Tool Comparison (AUROC)")

    # Panel B: Per-gene AUROC heatmap
    ax = axes[1]
    if tool_csv.exists():
        df = pd.read_csv(tool_csv)
        pivot = df.pivot_table(index="gene", columns="tool", values="auroc")
        if HAS_SNS:
            sns.heatmap(pivot, ax=ax, annot=True, fmt=".2f",
                       cmap="RdYlGn", vmin=0.5, vmax=1.0)
    else:
        ax.text(0.5, 0.5, "Data not yet\navailable",
                ha="center", va="center", transform=ax.transAxes)
    ax.set_title("B. Per-Gene AUROC")

    plt.tight_layout()
    out = FIG_DIR / "fig4_clinvar_validation.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    log.info("Saved: {}".format(out))


# =============================================================================
# Figure 6: Non-Coding Validation
# =============================================================================

def fig6_noncoding():
    """TERT MPRA scatter + ENCODE enrichment bars."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Evo2 vs TERT MPRA
    ax = axes[0]
    mpra_csv = RESULTS_DIR / "noncoding" / "tert_mpra_scatter.csv"
    if mpra_csv.exists():
        df = pd.read_csv(mpra_csv)
        ax.scatter(df["mpra_score"], df["evo2_score"],
                  alpha=0.3, s=8, color=COLORS["Evo2"])
        ax.set_xlabel("MPRA Expression Effect")
        ax.set_ylabel("Evo2 Score (-delta)")

        if len(df) > 10:
            from scipy import stats
            rho, p = stats.spearmanr(df["mpra_score"], df["evo2_score"])
            ax.text(0.05, 0.95, "rho = {:.3f}".format(rho),
                   transform=ax.transAxes, va="top", fontsize=10)
    else:
        ax.text(0.5, 0.5, "MPRA data not yet\navailable",
                ha="center", va="center", transform=ax.transAxes)
    ax.set_title("A. Evo2 vs TERT MPRA (Kircher 2019)")

    # Panel B: ENCODE enrichment
    ax = axes[1]
    encode_csv = RESULTS_DIR / "noncoding" / "encode_enrichment.csv"
    if encode_csv.exists():
        df = pd.read_csv(encode_csv)
        # Average across genes
        mean_enrich = df.groupby("ccre_class")["fold_enrichment"].mean()
        colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
        mean_enrich.plot(kind="bar", ax=ax, color=colors[:len(mean_enrich)])
        ax.axhline(y=1, color="k", linestyle="--", alpha=0.5)
        ax.set_ylabel("Fold Enrichment")
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=45)
    else:
        ax.text(0.5, 0.5, "ENCODE data not yet\navailable",
                ha="center", va="center", transform=ax.transAxes)
    ax.set_title("B. ENCODE cCRE Enrichment")

    plt.tight_layout()
    out = FIG_DIR / "fig6_noncoding.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    log.info("Saved: {}".format(out))


# =============================================================================
# Figure 7: Radiation Signature
# =============================================================================

def fig7_radiation():
    """Radiation-type vs other mutation impact."""
    json_file = RESULTS_DIR / "radiation_signatures" / "radiation_results.json"
    if not Path(json_file).exists():
        log.warning("No radiation data: {}".format(json_file))
        return

    with open(json_file) as f:
        results = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Per-gene radiation vs other mean delta
    ax = axes[0]
    genes = []
    rad_means = []
    other_means = []
    for r in results:
        comp = r.get("radiation_vs_other", {})
        if comp:
            genes.append(r["gene"])
            rad_means.append(comp["radiation_mean_delta"])
            other_means.append(comp["other_mean_delta"])

    if genes:
        x = np.arange(len(genes))
        width = 0.35
        ax.bar(x - width / 2, rad_means, width, label="Radiation (C>A/G>T)",
               color=COLORS["radiation"])
        ax.bar(x + width / 2, other_means, width, label="Other mutations",
               color=COLORS["other"])
        ax.set_xticks(x)
        ax.set_xticklabels(genes, rotation=45)
        ax.set_ylabel("Mean Delta")
        ax.legend()
    ax.set_title("A. Radiation vs Other Mutation Impact")

    # Panel B: Summary table of mutation classes
    ax = axes[1]
    csv_file = RESULTS_DIR / "radiation_signatures" / "radiation_class_summary.csv"
    if Path(csv_file).exists():
        df = pd.read_csv(csv_file)
        # Group by mutation class
        class_means = df.groupby("mutation_class")["mean_delta"].agg(["mean", "std"])
        class_means = class_means.sort_values("mean")
        colors_list = [COLORS.get("radiation", "#e41a1c")
                       if "oxidative" in c else COLORS["other"]
                       for c in class_means.index]
        class_means["mean"].plot(kind="barh", ax=ax, xerr=class_means["std"],
                                 color=colors_list, capsize=3)
        ax.set_xlabel("Mean Delta")
    else:
        ax.text(0.5, 0.5, "Data not yet\navailable",
                ha="center", va="center", transform=ax.transAxes)
    ax.set_title("B. Mutation Class Impact")

    plt.tight_layout()
    out = FIG_DIR / "fig7_radiation.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    log.info("Saved: {}".format(out))


# =============================================================================
# Figure 8: Astronaut Variants
# =============================================================================

def fig8_astronaut():
    """Astronaut variant scoring and context."""
    csv_file = RESULTS_DIR / "astronaut_variants" / "astronaut_variants.csv"
    if not Path(csv_file).exists():
        log.warning("No astronaut data: {}".format(csv_file))
        return

    df = pd.read_csv(csv_file)
    scored = df.dropna(subset=["evo2_delta"])

    if scored.empty:
        log.warning("No scored astronaut variants")
        return

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    # Lollipop plot of astronaut variant scores
    y_pos = np.arange(len(scored))
    colors = [COLORS.get("P/LP") if d < 0 else COLORS.get("B/LB")
              for d in scored["evo2_delta"]]

    ax.barh(y_pos, scored["evo2_delta"], color=colors, height=0.6)
    ax.set_yticks(y_pos)
    labels = ["{} {}".format(row["gene"], row["note"][:25])
              for _, row in scored.iterrows()]
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Evo2 Delta (more negative = more damaging)")
    ax.axvline(x=0, color="k", linestyle="-", alpha=0.3)
    ax.set_title("Astronaut Variant Evo2 Scores")

    plt.tight_layout()
    out = FIG_DIR / "fig8_astronaut.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    log.info("Saved: {}".format(out))


# =============================================================================
# Supplementary Figure S1: Orthogonality Matrix
# =============================================================================

def supp_s1_orthogonality():
    """Pairwise Spearman correlations between tools."""
    json_file = RESULTS_DIR / "benchmark_comparison" / "benchmark_results.json"
    if not Path(json_file).exists():
        return

    with open(json_file) as f:
        results = json.load(f)

    # Collect all pairwise correlations
    all_pairs = {}
    for r in results:
        for pair, data in r.get("orthogonality", {}).items():
            if pair not in all_pairs:
                all_pairs[pair] = []
            all_pairs[pair].append(data["rho"])

    if not all_pairs:
        return

    # Build matrix
    tools = set()
    for pair in all_pairs:
        t1, t2 = pair.split("_vs_")
        tools.add(t1)
        tools.add(t2)
    tools = sorted(tools)

    n = len(tools)
    matrix = np.eye(n)
    for pair, rhos in all_pairs.items():
        t1, t2 = pair.split("_vs_")
        i = tools.index(t1)
        j = tools.index(t2)
        mean_rho = np.mean(rhos)
        matrix[i, j] = mean_rho
        matrix[j, i] = mean_rho

    fig, ax = plt.subplots(figsize=(8, 6))
    if HAS_SNS:
        sns.heatmap(pd.DataFrame(matrix, index=tools, columns=tools),
                   ax=ax, annot=True, fmt=".2f", cmap="RdBu_r",
                   vmin=-1, vmax=1, center=0)
    ax.set_title("Tool Orthogonality (Mean Spearman rho)")

    plt.tight_layout()
    out = FIG_DIR / "supp_s1_orthogonality.pdf"
    fig.savefig(out)
    plt.close(fig)
    log.info("Saved: {}".format(out))


# =============================================================================
# Supplementary Figure S3: Indel Analysis
# =============================================================================

def supp_s3_indels():
    """In-frame vs frameshift indel analysis."""
    csv_file = RESULTS_DIR / "indel_analysis" / "indel_deltas.csv"
    if not Path(csv_file).exists():
        return

    df = pd.read_csv(csv_file)

    fig, ax = plt.subplots(figsize=(8, 5))

    if HAS_SNS:
        sns.boxplot(data=df, x="gene", y="mean_delta", hue="category",
                   ax=ax, palette={"frameshift": "#d62728", "inframe": "#1f77b4"})
    ax.set_ylabel("Mean Delta")
    ax.set_title("Frameshift vs In-Frame Indels")
    ax.legend(title="Category")

    plt.tight_layout()
    out = FIG_DIR / "supp_s3_indels.pdf"
    fig.savefig(out)
    plt.close(fig)
    log.info("Saved: {}".format(out))


# =============================================================================
# Main
# =============================================================================

FIGURE_MAP = {
    "2": ("Window Size Ablation", fig2_ablation),
    "3": ("DMS Calibration", fig3_dms_calibration),
    "4": ("ClinVar Validation", fig4_clinvar_validation),
    "6": ("Non-Coding Validation", fig6_noncoding),
    "7": ("Radiation Signatures", fig7_radiation),
    "8": ("Astronaut Variants", fig8_astronaut),
    "s1": ("Orthogonality Matrix", supp_s1_orthogonality),
    "s3": ("Indel Analysis", supp_s3_indels),
}


def main():
    parser = argparse.ArgumentParser(description="Generate publication figures")
    parser.add_argument("--fig", type=str, default="all",
                        help="Figure number (2-8, s1-s6) or 'all'")
    args = parser.parse_args()

    if not HAS_MPL:
        log.error("matplotlib required: pip install matplotlib seaborn")
        return

    setup_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    if args.fig == "all":
        figs_to_make = list(FIGURE_MAP.keys())
    else:
        figs_to_make = [args.fig]

    for fig_id in figs_to_make:
        if fig_id not in FIGURE_MAP:
            log.warning("Unknown figure: {}".format(fig_id))
            continue
        name, func = FIGURE_MAP[fig_id]
        log.info("Generating Fig {}: {}".format(fig_id, name))
        try:
            func()
        except Exception as e:
            log.error("Fig {} failed: {}".format(fig_id, e))

    log.info("\nFigures saved to {}".format(FIG_DIR))


if __name__ == "__main__":
    main()

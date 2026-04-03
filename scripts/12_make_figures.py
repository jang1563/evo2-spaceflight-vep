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
        json_file = RESULTS_DIR / gene / "calibration" / "dms_evo2_matched.json"
        cal_file = RESULTS_DIR / gene / "calibration" / "calibration_result.json"

        if json_file.exists():
            import json as _json
            matched = _json.load(open(json_file))
            df = pd.DataFrame(matched)
            ax.scatter(df["dms_score"], df["evo2_delta"],
                      alpha=0.3, s=8, color=COLORS["Evo2"])

            # Pull Spearman rho from calibration result if available
            rho_str = ""
            if cal_file.exists():
                cal = _json.load(open(cal_file))
                rho = cal.get("spearman_rho")
                p = cal.get("spearman_p", 1.0)
                if rho is not None:
                    rho_str = "rho = {:.3f}".format(rho)
                    if p < 0.001:
                        rho_str += "\np < 0.001"
            if not rho_str and len(df) > 10:
                from scipy import stats as _stats
                rho, p = _stats.spearmanr(df["dms_score"], df["evo2_delta"])
                rho_str = "rho = {:.3f}".format(rho)
            if rho_str:
                n_str = "n = {:,}".format(len(df))
                ax.text(0.05, 0.95, "{}\n{}".format(rho_str, n_str),
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
        pivot = df.pivot_table(index="gene", columns="tool", values="auroc",
                               aggfunc="first")
        if HAS_SNS:
            # Use custom annotation to show "—" for NaN instead of phantom values
            annot = pivot.copy()
            annot_text = annot.map(
                lambda x: "{:.3f}".format(x) if pd.notna(x) else "--")
            sns.heatmap(pivot, ax=ax, annot=annot_text, fmt="",
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
            p_str = "p = {:.1e}".format(p) if p < 0.001 else "p = {:.3f}".format(p)
            ax.text(0.05, 0.95,
                   "rho = {:.3f}\n{}\nn = {:,}".format(rho, p_str, len(df)),
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
        # Clean up raw variable names for display
        label_map = {
            "transversion_C_to_A": "C>A / G>T (oxidative)",
            "transversion_G_to_T": "C>A / G>T (oxidative)",
            "oxidative_C_to_A": "C>A / G>T (oxidative)",
            "radiation_oxidative": "C>A / G>T (oxidative)",
            "transition_C_to_T": "C>T / G>A (deamination)",
            "transition_cg_to_ta": "C>T / G>A (deamination)",
            "transition_CG_to_TA": "C>T / G>A (deamination)",
            "transition_AT_to_GC": "A>G / T>C (transition)",
            "transition_at_to_gc": "A>G / T>C (transition)",
            "transition_GC_to_AT": "G>A / C>T (deamination)",
            "transition_gc_to_at": "G>A / C>T (deamination)",
            "transversion_A_to_T": "A>T / T>A (transversion)",
            "transversion_A_to_C": "A>C / T>G (transversion)",
            "transversion_G_to_C": "G>C / C>G (transversion)",
            "other_transversion": "Other transversions",
        }
        clean_labels = [label_map.get(c, c.replace("_", " ").title())
                        for c in class_means.index]
        colors_list = [COLORS.get("radiation", "#e41a1c")
                       if "oxidative" in c or "C_to_A" in c or "G_to_T" in c
                       else COLORS["other"]
                       for c in class_means.index]
        class_means["mean"].plot(kind="barh", ax=ax, xerr=class_means["std"],
                                 color=colors_list, capsize=3)
        ax.set_yticklabels(clean_labels)
        ax.set_ylabel("")
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
    labels = ["{} {}".format(row["gene"], row.get("note", row.get("variant", "")))
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
# Figure 5: Per-Gene Constraint Landscapes
# =============================================================================

def _format_genomic_pos(x, _pos=None):
    """Format genomic position as Mb for tick labels."""
    return "{:.2f} Mb".format(x / 1e6)


def fig5_constraint_landscapes():
    """Per-gene constraint landscapes from Evo2 variant scores.

    For each gene, computes a positional constraint profile by averaging
    |delta| across all SNVs at each position, then smoothing with a
    rolling window. Exon structure is shown as a gene model track.
    """
    from matplotlib.ticker import FuncFormatter
    from utils.gene_coordinates import get_gene
    all_genes = [g.symbol for g in CONTROL_GENES] + [g.symbol for g in NOVEL_GENES]

    n_genes = len(all_genes)
    n_cols = 2
    n_rows = (n_genes + 1) // n_cols
    fig, axes_grid = plt.subplots(n_rows, n_cols, figsize=(14, 2.8 * n_rows),
                                   squeeze=False)
    axes = [axes_grid[i // n_cols, i % n_cols] for i in range(n_genes)]
    # Hide any unused axes
    for i in range(n_genes, n_rows * n_cols):
        axes_grid[i // n_cols, i % n_cols].set_visible(False)

    smooth_window = 50  # number of positions for rolling mean

    for idx, gene_symbol in enumerate(all_genes):
        ax = axes[idx]
        results_file = RESULTS_DIR / gene_symbol / "w8192" / "results.jsonl"

        if not results_file.exists():
            ax.text(0.5, 0.5, "{}: data not available".format(gene_symbol),
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_title(gene_symbol, fontsize=10, fontweight="bold")
            continue

        # Load SNV scores
        positions = []
        deltas = []
        clinvar_plp_pos = []
        with open(results_file) as f:
            for line in f:
                v = json.loads(line)
                # SNVs only (skip indels)
                if len(v.get("ref", "")) != 1 or len(v.get("alt", "")) != 1:
                    continue
                d = v.get("delta")
                if d is None:
                    continue
                positions.append(v["pos"])
                deltas.append(abs(d))
                # Track ClinVar P/LP positions
                cclass = v.get("clinvar_class", "")
                if cclass in ("P/LP", "Pathogenic", "Likely pathogenic",
                              "Pathogenic/Likely pathogenic"):
                    clinvar_plp_pos.append((v["pos"], abs(d)))

        if not positions:
            ax.text(0.5, 0.5, "{}: no SNV scores".format(gene_symbol),
                    ha="center", va="center", transform=ax.transAxes)
            continue

        # Build per-position mean |delta|
        pos_arr = np.array(positions)
        delta_arr = np.array(deltas)

        # Aggregate by position (mean across alt alleles at same position)
        pos_unique = np.unique(pos_arr)
        mean_delta = np.zeros(len(pos_unique))
        for i, p in enumerate(pos_unique):
            mask = pos_arr == p
            mean_delta[i] = delta_arr[mask].mean()

        # Rolling mean smoothing
        if len(mean_delta) > smooth_window:
            kernel = np.ones(smooth_window) / smooth_window
            smoothed = np.convolve(mean_delta, kernel, mode="same")
        else:
            smoothed = mean_delta

        # Plot constraint landscape
        ax.fill_between(pos_unique, 0, smoothed,
                        color=COLORS["Evo2"], alpha=0.4, linewidth=0)
        ax.plot(pos_unique, smoothed, color=COLORS["Evo2"],
                linewidth=0.5, alpha=0.8)

        # Mark ClinVar P/LP
        if clinvar_plp_pos:
            plp_x = [p[0] for p in clinvar_plp_pos]
            plp_y = [p[1] for p in clinvar_plp_pos]
            ax.scatter(plp_x, plp_y, color=COLORS["P/LP"],
                      s=6, alpha=0.5, zorder=3)

        # Gene structure annotation track
        try:
            gene_info = get_gene(gene_symbol)
            ymin, ymax = ax.get_ylim()
            track_y = ymin - (ymax - ymin) * 0.12
            track_h = (ymax - ymin) * 0.06

            # Gene body line
            ax.plot([gene_info.start, gene_info.end],
                    [track_y, track_y], color="black", linewidth=1.5,
                    clip_on=False)

            # Exons as filled rectangles
            for exon in gene_info.exons:
                is_coding = (exon.end > gene_info.cds_start and
                             exon.start < gene_info.cds_end)
                color = "#2c3e50" if is_coding else "#95a5a6"
                rect = plt.Rectangle(
                    (exon.start, track_y - track_h / 2),
                    exon.end - exon.start, track_h,
                    facecolor=color, edgecolor="black",
                    linewidth=0.5, clip_on=False, zorder=4)
                ax.add_patch(rect)
        except Exception:
            pass

        # Formatting
        ax.set_title(gene_symbol, fontsize=10, fontweight="bold")
        ax.tick_params(axis="y", labelsize=7)
        ax.tick_params(axis="x", labelsize=7)
        ax.xaxis.set_major_formatter(FuncFormatter(_format_genomic_pos))
        ax.set_ylabel("|delta|", fontsize=8)

        # Compact axis
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Add shared legend to first panel
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    legend_elements = [
        Patch(facecolor=COLORS["Evo2"], alpha=0.4,
              label="Evo2 |delta| (smoothed)"),
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=COLORS["P/LP"], markersize=5,
               label="ClinVar P/LP"),
        Patch(facecolor="#2c3e50", label="Coding exon"),
        Patch(facecolor="#95a5a6", label="Non-coding exon"),
    ]
    axes[0].legend(handles=legend_elements, loc="upper right",
                   fontsize=7, framealpha=0.8)

    fig.suptitle("Per-Gene Constraint Landscapes (Evo2)", fontsize=13, y=1.01)
    plt.tight_layout()
    out = FIG_DIR / "fig5_constraint_landscapes.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: {}".format(out))


# =============================================================================
# Main
# =============================================================================

FIGURE_MAP = {
    "2": ("Window Size Ablation", fig2_ablation),
    "3": ("DMS Calibration", fig3_dms_calibration),
    "4": ("ClinVar Validation", fig4_clinvar_validation),
    "5": ("Constraint Landscapes", fig5_constraint_landscapes),
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

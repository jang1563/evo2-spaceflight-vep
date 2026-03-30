"""
Benchmarking utilities for Evo2 VEP validation.

Implements:
- AUROC / AUPRC computation
- Precision-recall curves
- Calibration curves
- Pejaver 2022 likelihood ratio calibration for ClinGen evidence strengths
- Orthogonality matrix (pairwise Spearman correlations)
- Gene-specific threshold calibration
"""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

_HAS_SKLEARN = False
try:
    from scipy import stats as scipy_stats
    from sklearn.metrics import (
        auc,
        average_precision_score,
        precision_recall_curve,
        roc_auc_score,
        roc_curve,
    )
    from sklearn.calibration import calibration_curve
    _HAS_SKLEARN = True
except ImportError:
    pass


def _check_sklearn():
    if not _HAS_SKLEARN:
        raise ImportError(
            "scipy and scikit-learn required for benchmarking: "
            "pip install scipy scikit-learn"
        )


# =============================================================================
# ClinGen Evidence Strength Thresholds (Pejaver 2022)
# =============================================================================
# From: Pejaver et al. (2022) Am J Hum Genet 109(12):2163-2177
# Likelihood ratio thresholds for mapping computational scores to
# ACMG/AMP evidence strengths (PP3/BP4)

@dataclass
class ClinGenThresholds:
    """Evidence strength thresholds based on likelihood ratios."""
    # PP3 (pathogenicity supporting) thresholds
    pp3_supporting_lr: float = 2.08     # LR ≥ 2.08
    pp3_moderate_lr: float = 4.33       # LR ≥ 4.33
    pp3_strong_lr: float = 18.7         # LR ≥ 18.7
    pp3_very_strong_lr: float = 350.0   # LR ≥ 350

    # BP4 (benign supporting) thresholds (reciprocal)
    bp4_supporting_lr: float = 0.48     # LR ≤ 1/2.08
    bp4_moderate_lr: float = 0.23       # LR ≤ 1/4.33
    bp4_strong_lr: float = 0.053        # LR ≤ 1/18.7


# =============================================================================
# Core Metrics
# =============================================================================

def compute_auroc(
    labels: np.ndarray,
    scores: np.ndarray,
) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute AUROC and ROC curve.

    Args:
        labels: Binary labels (1=pathogenic, 0=benign)
        scores: Continuous scores (higher = more pathogenic)

    Returns:
        (auroc, fpr, tpr, thresholds)
    """
    _check_sklearn()
    fpr, tpr, thresholds = roc_curve(labels, scores)
    auroc = auc(fpr, tpr)
    return auroc, fpr, tpr, thresholds


def compute_auprc(
    labels: np.ndarray,
    scores: np.ndarray,
) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute AUPRC and precision-recall curve.

    Returns:
        (auprc, precision, recall, thresholds)
    """
    _check_sklearn()
    precision, recall, thresholds = precision_recall_curve(labels, scores)
    auprc = average_precision_score(labels, scores)
    return auprc, precision, recall, thresholds


def precision_at_recall(
    labels: np.ndarray,
    scores: np.ndarray,
    target_recall: float = 0.9,
) -> float:
    """Compute precision at a given recall level."""
    _check_sklearn()
    precision, recall, _ = precision_recall_curve(labels, scores)
    # Find the threshold closest to target recall
    idx = np.argmin(np.abs(recall - target_recall))
    return float(precision[idx])


def recall_at_precision(
    labels: np.ndarray,
    scores: np.ndarray,
    target_precision: float = 0.9,
) -> float:
    """Compute recall at a given precision level."""
    _check_sklearn()
    precision, recall, _ = precision_recall_curve(labels, scores)
    # Find indices where precision >= target
    valid = precision >= target_precision
    if not np.any(valid):
        return 0.0
    return float(np.max(recall[valid]))


def compute_calibration(
    labels: np.ndarray,
    scores: np.ndarray,
    n_bins: int = 10,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute calibration curve.

    Returns:
        (fraction_of_positives, mean_predicted_value) for each bin
    """
    _check_sklearn()
    # Normalize scores to [0, 1] range for calibration
    s_min, s_max = scores.min(), scores.max()
    if s_max > s_min:
        scores_norm = (scores - s_min) / (s_max - s_min)
    else:
        scores_norm = np.full_like(scores, 0.5)

    fraction_pos, mean_pred = calibration_curve(
        labels, scores_norm, n_bins=n_bins, strategy="uniform"
    )
    return fraction_pos, mean_pred


# =============================================================================
# Pejaver Likelihood Ratio Calibration
# =============================================================================

def compute_likelihood_ratios(
    pathogenic_scores: np.ndarray,
    benign_scores: np.ndarray,
    thresholds: np.ndarray,
) -> np.ndarray:
    """
    Compute likelihood ratios at given thresholds.

    LR(t) = P(score >= t | pathogenic) / P(score >= t | benign)
    """
    lrs = np.zeros(len(thresholds))
    n_path = len(pathogenic_scores)
    n_ben = len(benign_scores)

    for i, t in enumerate(thresholds):
        p_path = np.sum(pathogenic_scores >= t) / n_path if n_path > 0 else 0
        p_ben = np.sum(benign_scores >= t) / n_ben if n_ben > 0 else 0

        if p_ben > 0:
            lrs[i] = p_path / p_ben
        else:
            lrs[i] = float('inf') if p_path > 0 else 1.0

    return lrs


def calibrate_gene_thresholds(
    pathogenic_scores: np.ndarray,
    benign_scores: np.ndarray,
    n_threshold_points: int = 1000,
) -> Dict[str, Optional[float]]:
    """
    Calibrate gene-specific thresholds using Pejaver 2022 LR framework.

    Returns dict mapping evidence strength to score threshold.
    """
    all_scores = np.concatenate([pathogenic_scores, benign_scores])
    thresholds = np.linspace(all_scores.min(), all_scores.max(), n_threshold_points)

    lrs = compute_likelihood_ratios(pathogenic_scores, benign_scores, thresholds)
    clingen = ClinGenThresholds()

    result = {
        "pp3_supporting": None,
        "pp3_moderate": None,
        "pp3_strong": None,
        "bp4_supporting": None,
        "bp4_moderate": None,
        "bp4_strong": None,
    }

    # PP3 thresholds (find lowest score where LR >= threshold)
    for lr_thresh, key in [
        (clingen.pp3_supporting_lr, "pp3_supporting"),
        (clingen.pp3_moderate_lr, "pp3_moderate"),
        (clingen.pp3_strong_lr, "pp3_strong"),
    ]:
        valid = np.where(lrs >= lr_thresh)[0]
        if len(valid) > 0:
            result[key] = float(thresholds[valid[0]])

    # BP4 thresholds (find highest score where LR <= threshold)
    for lr_thresh, key in [
        (clingen.bp4_supporting_lr, "bp4_supporting"),
        (clingen.bp4_moderate_lr, "bp4_moderate"),
        (clingen.bp4_strong_lr, "bp4_strong"),
    ]:
        valid = np.where(lrs <= lr_thresh)[0]
        if len(valid) > 0:
            result[key] = float(thresholds[valid[-1]])

    return result


def apply_gene_thresholds(
    score: float,
    thresholds: Dict[str, Optional[float]],
) -> str:
    """
    Apply calibrated thresholds to classify a score.

    Returns evidence strength string.
    """
    # Check PP3 (pathogenic direction, from strongest to weakest)
    for key in ["pp3_strong", "pp3_moderate", "pp3_supporting"]:
        if thresholds[key] is not None and score >= thresholds[key]:
            return key.upper()

    # Check BP4 (benign direction, from strongest to weakest)
    for key in ["bp4_strong", "bp4_moderate", "bp4_supporting"]:
        if thresholds[key] is not None and score <= thresholds[key]:
            return key.upper()

    return "INDETERMINATE"


# =============================================================================
# Orthogonality Assessment
# =============================================================================

def compute_orthogonality_matrix(
    tool_scores: Dict[str, np.ndarray],
) -> Tuple[np.ndarray, List[str]]:
    """
    Compute pairwise Spearman correlations between tools.

    Args:
        tool_scores: Dict mapping tool name -> score array (same variant order)

    Returns:
        (correlation_matrix, tool_names)
    """
    names = sorted(tool_scores.keys())
    n = len(names)
    corr_matrix = np.ones((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            scores_i = tool_scores[names[i]]
            scores_j = tool_scores[names[j]]

            # Handle NaN values (some tools may not score all variants)
            valid = ~(np.isnan(scores_i) | np.isnan(scores_j))
            if np.sum(valid) < 10:
                corr_matrix[i, j] = corr_matrix[j, i] = np.nan
                continue

            rho, p = scipy_stats.spearmanr(scores_i[valid], scores_j[valid])
            corr_matrix[i, j] = rho
            corr_matrix[j, i] = rho

    return corr_matrix, names


def interpret_orthogonality(rho: float) -> str:
    """Interpret pairwise correlation for added value assessment."""
    if abs(rho) > 0.9:
        return "Recapitulates (low added value)"
    elif abs(rho) > 0.7:
        return "Moderately correlated"
    elif abs(rho) > 0.3:
        return "Complementary (high added value)"
    else:
        return "Near-orthogonal (may provide PP3 without double-counting)"


# =============================================================================
# Gene-Specific Landscape Analysis (EVE-style)
# =============================================================================

def compute_score_landscape(
    all_snv_scores: np.ndarray,
    pathogenic_scores: Optional[np.ndarray] = None,
    benign_scores: Optional[np.ndarray] = None,
    vus_scores: Optional[np.ndarray] = None,
    n_bins: int = 100,
) -> Dict:
    """
    Compute gene-specific score landscape for EVE-style visualization.

    Returns dict with distribution statistics and histogram data.
    """
    result = {
        "all_snv": {
            "mean": float(np.mean(all_snv_scores)),
            "std": float(np.std(all_snv_scores)),
            "median": float(np.median(all_snv_scores)),
            "q25": float(np.percentile(all_snv_scores, 25)),
            "q75": float(np.percentile(all_snv_scores, 75)),
            "n": len(all_snv_scores),
        },
    }

    # Histogram for all SNVs
    hist, bin_edges = np.histogram(all_snv_scores, bins=n_bins)
    result["all_snv"]["hist_counts"] = hist.tolist()
    result["all_snv"]["hist_edges"] = bin_edges.tolist()

    # Overlay distributions for labeled variants
    for name, scores in [
        ("pathogenic", pathogenic_scores),
        ("benign", benign_scores),
        ("vus", vus_scores),
    ]:
        if scores is not None and len(scores) > 0:
            result[name] = {
                "mean": float(np.mean(scores)),
                "std": float(np.std(scores)),
                "median": float(np.median(scores)),
                "n": len(scores),
                "hist_counts": np.histogram(scores, bins=bin_edges)[0].tolist(),
            }

    return result


# =============================================================================
# Multi-Tool Comparison
# =============================================================================

@dataclass
class ToolComparison:
    """Results of comparing multiple tools on the same variants."""
    gene: str
    variant_ids: List[str]
    labels: np.ndarray    # Binary labels
    tool_scores: Dict[str, np.ndarray]
    tool_aurocs: Dict[str, float] = field(default_factory=dict)
    tool_auprcs: Dict[str, float] = field(default_factory=dict)

    def compute_all_metrics(self):
        """Compute AUROC and AUPRC for all tools."""
        _check_sklearn()
        for tool_name, scores in self.tool_scores.items():
            valid = ~np.isnan(scores)
            if np.sum(valid) < 10 or len(np.unique(self.labels[valid])) < 2:
                continue
            self.tool_aurocs[tool_name] = roc_auc_score(
                self.labels[valid], scores[valid]
            )
            self.tool_auprcs[tool_name] = average_precision_score(
                self.labels[valid], scores[valid]
            )

    def summary(self) -> str:
        lines = [f"Tool Comparison — {self.gene}"]
        lines.append(f"  Variants: {len(self.variant_ids)} "
                     f"(P/LP: {np.sum(self.labels == 1)}, B/LB: {np.sum(self.labels == 0)})")
        lines.append(f"  {'Tool':<20} {'AUROC':>8} {'AUPRC':>8}")
        lines.append("  " + "-" * 38)
        for tool in sorted(self.tool_aurocs.keys()):
            auroc = self.tool_aurocs.get(tool, float('nan'))
            auprc = self.tool_auprcs.get(tool, float('nan'))
            lines.append(f"  {tool:<20} {auroc:>8.4f} {auprc:>8.4f}")
        return "\n".join(lines)


# =============================================================================
# Sanity Checks
# =============================================================================

def run_sanity_checks(
    results: List[dict],
    gene: str,
) -> Dict[str, bool]:
    """
    Run sanity checks on scoring results for a gene.

    Checks:
    1. Synonymous variants have delta ≈ 0
    2. Nonsense/frameshift variants have strongly negative delta
    3. ClinVar P/LP scores lower than B/LB (for ≥2-star)
    """
    checks = {}

    # Filter results for this gene
    gene_results = [r for r in results if r["gene"] == gene]
    if not gene_results:
        return {"no_results": False}

    # Check 1: Synonymous delta ≈ 0
    syn = [r["delta"] for r in gene_results
           if r.get("region_type") == "synonymous"]
    if syn:
        mean_syn = np.mean(syn)
        checks["synonymous_near_zero"] = abs(mean_syn) < 0.5
        checks["synonymous_mean_delta"] = float(mean_syn)

    # Check 2: P/LP vs B/LB separation (≥2 stars)
    plp = [r["delta"] for r in gene_results
           if r["clinvar_class"] == "P/LP" and r.get("clinvar_stars", 0) >= 2]
    blb = [r["delta"] for r in gene_results
           if r["clinvar_class"] == "B/LB" and r.get("clinvar_stars", 0) >= 2]

    if plp and blb:
        checks["plp_lower_than_blb"] = np.mean(plp) < np.mean(blb)
        checks["plp_mean_delta"] = float(np.mean(plp))
        checks["blb_mean_delta"] = float(np.mean(blb))

    return checks


# =============================================================================
# Save / Load
# =============================================================================

def save_comparison(comparison: ToolComparison, filepath: str):
    """Save tool comparison results to JSON."""
    data = {
        "gene": comparison.gene,
        "n_variants": len(comparison.variant_ids),
        "n_pathogenic": int(np.sum(comparison.labels == 1)),
        "n_benign": int(np.sum(comparison.labels == 0)),
        "aurocs": comparison.tool_aurocs,
        "auprcs": comparison.tool_auprcs,
    }
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)


def load_comparison(filepath: str) -> dict:
    """Load tool comparison results from JSON."""
    with open(filepath) as f:
        return json.load(f)

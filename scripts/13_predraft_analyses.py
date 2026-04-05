#!/usr/bin/env python3
"""
Pre-draft manuscript analyses (Phase 8b).

Runs ALL analyses required before manuscript drafting:
  1. Matched-variant benchmarking (intersection of Evo2 ∩ CADD ∩ AM ∩ REVEL)
  2. Bootstrap 95% CIs on all AUROCs
  3. DeLong's test for pairwise AUROC comparisons
  4. DMS multi-tool comparison (CADD/AM/REVEL vs DMS on same variants)
  5. Transversion control (radiation-oxidative vs other transversions)
  6. Ensemble model (Evo2 + CADD logistic regression, 5-fold CV)
  7. MPRA multi-tool comparison (TERT promoter)
  8. Regenerate missing calibration files (BRCA1, CHEK2)
  9. Regenerate indel AUROCs for all genes
 10. Exact total variant count

Usage:
  python 13_predraft_analyses.py --analysis all
  python 13_predraft_analyses.py --analysis matched_benchmark
  python 13_predraft_analyses.py --analysis dms_multitool
  python 13_predraft_analyses.py --analysis transversion_control
  python 13_predraft_analyses.py --analysis ensemble
  python 13_predraft_analyses.py --analysis mpra_multitool
  python 13_predraft_analyses.py --analysis regenerate_calibration
  python 13_predraft_analyses.py --analysis regenerate_indels
  python 13_predraft_analyses.py --analysis variant_count
"""

import argparse
import csv
import gzip
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from utils.config import (
    ALPHAMISSENSE_TSV,
    CLINVAR_VCF,
    DMS_DIR,
    MPRA_DIR,
    RESULTS_DIR,
    SHARED_DATA_DIR,
)
from utils.gene_coordinates import CONTROL_GENES, GENES, get_gene
from utils.clinvar_parser import parse_clinvar_vcf, clinvar_to_variants
from utils.sequence_utils import ScoringCheckpoint
from utils.benchmarking import (
    compute_auroc,
    compute_auprc,
    interpret_orthogonality,
)
from utils.hgvs_mapping import (
    map_brca1_dms,
    map_chek2_dms,
    map_dnmt3a_garcia_dms,
    map_tp53_dms,
    map_tert_mpra,
)

from scipy import stats as scipy_stats
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

OUT_DIR = RESULTS_DIR / "predraft"

# Genes with sufficient ClinVar variants for benchmarking
BENCHMARK_GENES = ["ATM", "BRCA1", "TP53", "CHEK2", "DNMT3A", "TERT"]
# Genes with DMS data
DMS_GENES = ["BRCA1", "TP53", "CHEK2", "DNMT3A"]
# Tool comparison set for matched benchmarking
MATCHED_TOOLS = ["Evo2", "CADD", "AlphaMissense", "REVEL"]


# =============================================================================
# Helper: Load tool scores (reused from 07_benchmark_tools.py)
# =============================================================================

def load_tool_scores_for_gene(gene_symbol: str, variant_keys: set) -> Dict[str, Dict[str, float]]:
    """Load scores from all benchmark tools for a set of variants."""
    gene = get_gene(gene_symbol)
    tool_scores = {}

    import pysam

    # AlphaMissense
    am_path = str(ALPHAMISSENSE_TSV)
    if Path(am_path).exists():
        scores = {}
        tbx = pysam.TabixFile(am_path)
        chrom_q = gene.chrom if gene.chrom.startswith("chr") else "chr" + gene.chrom
        try:
            for row in tbx.fetch(chrom_q, gene.start, gene.end):
                fields = row.split("\t")
                chrom = fields[0] if fields[0].startswith("chr") else "chr" + fields[0]
                pos = int(fields[1]) - 1
                key = f"{chrom}:{pos}:{fields[2]}>{fields[3]}"
                if key in variant_keys:
                    scores[key] = float(fields[8])
        except Exception:
            pass
        tbx.close()
        if scores:
            tool_scores["AlphaMissense"] = scores

    # CADD v1.7
    cadd_path = str(SHARED_DATA_DIR / "benchmarks" / "cadd" / "whole_genome_SNVs.tsv.gz")
    if Path(cadd_path).exists():
        scores = {}
        tbx = pysam.TabixFile(cadd_path)
        chrom_q = gene.chrom.replace("chr", "")
        try:
            for row in tbx.fetch(chrom_q, gene.start, gene.end):
                fields = row.split("\t")
                chrom = "chr" + fields[0]
                pos = int(fields[1]) - 1
                key = f"{chrom}:{pos}:{fields[2]}>{fields[3]}"
                if key in variant_keys:
                    scores[key] = float(fields[5])  # PHRED score
        except Exception:
            pass
        tbx.close()
        if scores:
            tool_scores["CADD"] = scores

    # REVEL
    revel_path = str(SHARED_DATA_DIR / "benchmarks" / "revel" / "revel_genes.bed.gz")
    if Path(revel_path).exists():
        scores = {}
        tbx = pysam.TabixFile(revel_path)
        chrom_q = gene.chrom if gene.chrom.startswith("chr") else "chr" + gene.chrom
        try:
            for row in tbx.fetch(chrom_q, gene.start, gene.end):
                fields = row.split("\t")
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                key = f"{chrom}:{pos}:{ref}>{alt}"
                if key in variant_keys:
                    scores[key] = float(fields[7])
        except Exception:
            pass
        tbx.close()
        if scores:
            tool_scores["REVEL"] = scores

    # SpliceAI
    spliceai_path = SHARED_DATA_DIR / "benchmarks" / "spliceai" / f"{gene_symbol}_spliceai_scores.tsv.gz"
    if spliceai_path.exists():
        scores = {}
        with gzip.open(str(spliceai_path), "rt") as fh:
            fh.readline()  # skip header
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                chrom = fields[0]
                pos = int(fields[1]) - 1
                key = f"{chrom}:{pos}:{fields[2]}>{fields[3]}"
                if key in variant_keys:
                    scores[key] = max(float(fields[5]), float(fields[6]),
                                      float(fields[7]), float(fields[8]))
        if scores:
            tool_scores["SpliceAI"] = scores

    # ncER
    ncer_path = SHARED_DATA_DIR / "noncoding" / "ncer" / "ncER_hg38_gene_regions.bed.gz"
    if ncer_path.exists():
        scores = {}
        tbx = pysam.TabixFile(str(ncer_path))
        chrom_q = gene.chrom if gene.chrom.startswith("chr") else "chr" + gene.chrom
        bin_scores = {}
        try:
            for row in tbx.fetch(chrom_q, gene.start, gene.end):
                fields = row.split("\t")
                bin_start = int(fields[1])
                bin_end = int(fields[2])
                ncer_score = float(fields[3])
                for p in range(bin_start, bin_end):
                    bin_scores[p] = ncer_score
        except Exception:
            pass
        tbx.close()
        for key in variant_keys:
            parts = key.split(":")
            pos = int(parts[1])
            if pos in bin_scores:
                scores[key] = bin_scores[pos]
        if scores:
            tool_scores["ncER"] = scores

    # LINSIGHT (via pyBigWig, hg19 lifted) — skip if chain file missing
    try:
        import pyBigWig
        from pyliftover import LiftOver
        linsight_path = SHARED_DATA_DIR / "benchmarks" / "linsight" / "LINSIGHT.bw"
        chain_path = SHARED_DATA_DIR / "liftover" / "hg38ToHg19.over.chain.gz"
        if linsight_path.exists() and chain_path.exists():
            lo = LiftOver(str(chain_path))
            bw = pyBigWig.open(str(linsight_path))
            scores = {}
            for key in variant_keys:
                parts = key.split(":")
                chrom38 = parts[0]
                pos38 = int(parts[1])
                result = lo.convert_coordinate(chrom38, pos38)
                if result and len(result) > 0:
                    chrom19, pos19 = result[0][0], result[0][1]
                    try:
                        val = bw.values(chrom19, pos19, pos19 + 1)
                        if val and val[0] is not None:
                            scores[key] = val[0]
                    except Exception:
                        pass
            bw.close()
            if scores:
                tool_scores["LINSIGHT"] = scores
    except (ImportError, FileNotFoundError):
        pass

    return tool_scores


def load_evo2_deltas(gene_symbol: str, window_size: int = 8192) -> Dict[str, float]:
    """Load Evo2 delta scores for a gene."""
    ckpt_dir = RESULTS_DIR / gene_symbol / f"w{window_size}"
    checkpoint = ScoringCheckpoint(str(ckpt_dir))
    results = checkpoint.load_results()
    deltas = {}
    for r in results:
        key = f"{r['chrom']}:{r['pos']}:{r['ref']}>{r['alt']}"
        deltas[key] = r["delta"]
    return deltas


def get_clinvar_variants(gene_symbol: str, min_stars: int = 2):
    """Get ClinVar P/LP and B/LB variant keys for a gene."""
    entries_plp = parse_clinvar_vcf(
        str(CLINVAR_VCF), genes={gene_symbol},
        min_stars=min_stars, significance_filter={"P/LP"},
    )
    entries_blb = parse_clinvar_vcf(
        str(CLINVAR_VCF), genes={gene_symbol},
        min_stars=min_stars, significance_filter={"B/LB"},
    )
    plp_keys = []
    blb_keys = []
    for e in clinvar_to_variants(entries_plp):
        if e.is_snv:
            plp_keys.append(e.key())
    for e in clinvar_to_variants(entries_blb):
        if e.is_snv:
            blb_keys.append(e.key())
    return plp_keys, blb_keys


# =============================================================================
# Helper: Bootstrap CI
# =============================================================================

def bootstrap_auroc(labels, scores, n_bootstrap=1000, ci=0.95, seed=42):
    """Compute AUROC with bootstrap 95% CI."""
    rng = np.random.RandomState(seed)
    n = len(labels)
    labels = np.asarray(labels)
    scores = np.asarray(scores)

    # Point estimate
    point = roc_auc_score(labels, scores)

    # Stratified bootstrap
    pos_idx = np.where(labels == 1)[0]
    neg_idx = np.where(labels == 0)[0]

    boot_aurocs = []
    for _ in range(n_bootstrap):
        p_sample = rng.choice(pos_idx, size=len(pos_idx), replace=True)
        n_sample = rng.choice(neg_idx, size=len(neg_idx), replace=True)
        idx = np.concatenate([p_sample, n_sample])
        try:
            ba = roc_auc_score(labels[idx], scores[idx])
            boot_aurocs.append(ba)
        except ValueError:
            continue

    boot_aurocs = np.array(boot_aurocs)
    alpha = (1 - ci) / 2
    ci_low = np.percentile(boot_aurocs, 100 * alpha)
    ci_high = np.percentile(boot_aurocs, 100 * (1 - alpha))

    return point, ci_low, ci_high, boot_aurocs


# =============================================================================
# Helper: DeLong test
# =============================================================================

def delong_test(labels, scores1, scores2):
    """
    DeLong's test for comparing two AUROCs on the same sample.
    Returns z-statistic and two-sided p-value.
    Based on Sun & Xu (2014) fast algorithm.
    """
    labels = np.asarray(labels)
    scores1 = np.asarray(scores1)
    scores2 = np.asarray(scores2)

    pos = np.where(labels == 1)[0]
    neg = np.where(labels == 0)[0]
    m = len(pos)
    n = len(neg)

    # Structural components for both ROC curves
    def compute_placements(scores):
        V10 = np.zeros(m)
        V01 = np.zeros(n)
        for i, pi in enumerate(pos):
            V10[i] = np.mean(scores[pi] > scores[neg]) + 0.5 * np.mean(scores[pi] == scores[neg])
        for j, nj in enumerate(neg):
            V01[j] = np.mean(scores[pos] > scores[nj]) + 0.5 * np.mean(scores[pos] == scores[nj])
        return V10, V01

    V10_1, V01_1 = compute_placements(scores1)
    V10_2, V01_2 = compute_placements(scores2)

    auc1 = np.mean(V10_1)
    auc2 = np.mean(V10_2)

    # Covariance matrix of (AUC1, AUC2)
    S10 = np.cov(np.stack([V10_1, V10_2])) / m if m > 1 else np.zeros((2, 2))
    S01 = np.cov(np.stack([V01_1, V01_2])) / n if n > 1 else np.zeros((2, 2))
    S = S10 + S01

    # z-test
    diff = auc1 - auc2
    var_diff = S[0, 0] + S[1, 1] - 2 * S[0, 1]
    if var_diff <= 0:
        return 0.0, 1.0
    z = diff / np.sqrt(var_diff)
    p = 2 * scipy_stats.norm.sf(abs(z))

    return float(z), float(p)


# =============================================================================
# Analysis 1: Matched-variant benchmarking + Bootstrap CIs + DeLong
# =============================================================================

def run_matched_benchmark():
    """
    Compute AUROCs on intersection of variants scored by ALL tools.
    Also compute per-tool AUROCs with bootstrap CIs and DeLong tests.
    """
    log.info("=" * 80)
    log.info("ANALYSIS 1-3: Matched-variant benchmarking + CIs + DeLong")
    log.info("=" * 80)

    all_results = []

    for gene_symbol in BENCHMARK_GENES:
        log.info(f"\n--- {gene_symbol} ---")

        # Load Evo2 scores
        evo2_deltas = load_evo2_deltas(gene_symbol)
        if not evo2_deltas:
            log.warning(f"  No Evo2 results for {gene_symbol}")
            continue

        # Get ClinVar variants
        plp_keys, blb_keys = get_clinvar_variants(gene_symbol)
        all_keys = plp_keys + blb_keys
        all_labels = [1] * len(plp_keys) + [0] * len(blb_keys)

        # Filter to Evo2-scored
        scored_mask = [k in evo2_deltas for k in all_keys]
        all_keys = [k for k, m in zip(all_keys, scored_mask) if m]
        all_labels = [l for l, m in zip(all_labels, scored_mask) if m]

        if sum(all_labels) < 5 or sum(1 - np.array(all_labels)) < 5:
            log.warning(f"  Too few ClinVar variants")
            continue

        # Load all tool scores
        variant_key_set = set(all_keys)
        tool_scores = load_tool_scores_for_gene(gene_symbol, variant_key_set)

        # Add Evo2 (negated: higher = more pathogenic)
        tool_scores["Evo2"] = {k: -evo2_deltas[k] for k in all_keys if k in evo2_deltas}

        log.info(f"  ClinVar ≥2★: {sum(all_labels)} P/LP, {len(all_labels) - sum(all_labels)} B/LB")
        for tn, ts in tool_scores.items():
            n_scored = sum(1 for k in all_keys if k in ts)
            log.info(f"  {tn}: {n_scored} variants scored")

        # --- Per-tool AUROCs with bootstrap CIs ---
        per_tool_results = {}
        per_tool_arrays = {}  # for DeLong

        for tool_name, scores_dict in tool_scores.items():
            t_labels = []
            t_scores = []
            for key, label in zip(all_keys, all_labels):
                if key in scores_dict:
                    t_labels.append(label)
                    t_scores.append(scores_dict[key])

            if len(t_labels) < 10 or len(set(t_labels)) < 2:
                continue

            t_labels = np.array(t_labels)
            t_scores = np.array(t_scores)

            auroc, ci_low, ci_high, _ = bootstrap_auroc(t_labels, t_scores)
            auprc_val, _, _, _ = compute_auprc(t_labels, t_scores)

            per_tool_results[tool_name] = {
                "auroc": float(auroc),
                "auroc_ci_low": float(ci_low),
                "auroc_ci_high": float(ci_high),
                "auprc": float(auprc_val),
                "n_scored": len(t_labels),
                "n_plp": int(sum(t_labels)),
                "n_blb": int(sum(1 - t_labels)),
            }
            per_tool_arrays[tool_name] = (t_labels, t_scores)

            log.info(f"  {tool_name}: AUROC={auroc:.4f} [{ci_low:.4f}, {ci_high:.4f}] (n={len(t_labels)})")

        # --- Matched-variant intersection (Evo2 ∩ CADD ∩ AM ∩ REVEL) ---
        matched_keys = []
        matched_labels = []
        available_tools = [t for t in MATCHED_TOOLS if t in tool_scores]

        for key, label in zip(all_keys, all_labels):
            if all(key in tool_scores[t] for t in available_tools):
                matched_keys.append(key)
                matched_labels.append(label)

        matched_results = {}
        matched_arrays = {}

        if len(matched_keys) >= 10 and len(set(matched_labels)) == 2:
            matched_labels_arr = np.array(matched_labels)
            log.info(f"\n  MATCHED SET: {len(matched_keys)} variants ({sum(matched_labels)} P/LP, {len(matched_labels) - sum(matched_labels)} B/LB)")
            log.info(f"  Tools in matched set: {available_tools}")

            for tool_name in available_tools:
                t_scores = np.array([tool_scores[tool_name][k] for k in matched_keys])
                auroc, ci_low, ci_high, _ = bootstrap_auroc(matched_labels_arr, t_scores)
                auprc_val, _, _, _ = compute_auprc(matched_labels_arr, t_scores)

                matched_results[tool_name] = {
                    "auroc": float(auroc),
                    "auroc_ci_low": float(ci_low),
                    "auroc_ci_high": float(ci_high),
                    "auprc": float(auprc_val),
                    "n": len(matched_keys),
                }
                matched_arrays[tool_name] = (matched_labels_arr, t_scores)
                log.info(f"  {tool_name} (matched): AUROC={auroc:.4f} [{ci_low:.4f}, {ci_high:.4f}]")

        # --- DeLong pairwise tests on matched set ---
        delong_results = {}
        if matched_arrays:
            tool_list = sorted(matched_arrays.keys())
            for i, t1 in enumerate(tool_list):
                for j, t2 in enumerate(tool_list):
                    if j <= i:
                        continue
                    labels1, scores1 = matched_arrays[t1]
                    _, scores2 = matched_arrays[t2]
                    z, p = delong_test(labels1, scores1, scores2)
                    pair = f"{t1}_vs_{t2}"
                    delong_results[pair] = {
                        "z_statistic": float(z),
                        "p_value": float(p),
                        "auroc_1": matched_results[t1]["auroc"],
                        "auroc_2": matched_results[t2]["auroc"],
                    }
                    log.info(f"  DeLong {t1} vs {t2}: z={z:.3f}, p={p:.4e}")

        # --- Pairwise orthogonality ---
        ortho = {}
        tool_list = sorted(tool_scores.keys())
        for i, t1 in enumerate(tool_list):
            for j, t2 in enumerate(tool_list):
                if j <= i:
                    continue
                shared = [k for k in all_keys if k in tool_scores[t1] and k in tool_scores[t2]]
                if len(shared) < 10:
                    continue
                s1 = np.array([tool_scores[t1][k] for k in shared])
                s2 = np.array([tool_scores[t2][k] for k in shared])
                rho, p = scipy_stats.spearmanr(s1, s2)
                ortho[f"{t1}_vs_{t2}"] = {
                    "spearman_rho": float(rho),
                    "p_value": float(p),
                    "n_shared": len(shared),
                }

        all_results.append({
            "gene": gene_symbol,
            "per_tool": per_tool_results,
            "matched_set": matched_results,
            "matched_n": len(matched_keys),
            "matched_n_plp": sum(matched_labels) if matched_labels else 0,
            "matched_n_blb": (len(matched_labels) - sum(matched_labels)) if matched_labels else 0,
            "delong": delong_results,
            "orthogonality": ortho,
        })

    # Save
    out_dir = OUT_DIR / "matched_benchmark"
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(out_dir / "matched_benchmark_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    # Save flat CSV for plotting
    rows = []
    for r in all_results:
        for tool, m in r["per_tool"].items():
            rows.append({
                "gene": r["gene"], "tool": tool, "set": "per_tool",
                "auroc": m["auroc"], "ci_low": m["auroc_ci_low"],
                "ci_high": m["auroc_ci_high"], "n": m["n_scored"],
            })
        for tool, m in r["matched_set"].items():
            rows.append({
                "gene": r["gene"], "tool": tool, "set": "matched",
                "auroc": m["auroc"], "ci_low": m["auroc_ci_low"],
                "ci_high": m["auroc_ci_high"], "n": m["n"],
            })
    pd.DataFrame(rows).to_csv(out_dir / "all_aurocs_with_ci.csv", index=False)

    log.info(f"\nSaved to {out_dir}")
    return all_results


# =============================================================================
# Analysis 4: DMS multi-tool comparison
# =============================================================================

def run_dms_multitool():
    """Compute DMS Spearman rho for all tools on same matched variants."""
    log.info("=" * 80)
    log.info("ANALYSIS 4: DMS multi-tool comparison")
    log.info("=" * 80)

    DMS_LOADERS = {
        "BRCA1": ("brca1", map_brca1_dms),
        "TP53": ("tp53", map_tp53_dms),
        "CHEK2": ("chek2", map_chek2_dms),
        "DNMT3A": ("dnmt3a", map_dnmt3a_garcia_dms),
    }

    all_results = []

    for gene_symbol in DMS_GENES:
        log.info(f"\n--- {gene_symbol} ---")

        # Load Evo2 scores
        evo2_deltas = load_evo2_deltas(gene_symbol)
        if not evo2_deltas:
            continue

        # Load DMS data
        subdir, mapper_func = DMS_LOADERS[gene_symbol]
        dms_dir = DMS_DIR / subdir

        # DMS file paths
        DMS_FILES = {
            "BRCA1": "brca1_findlay_sge.csv",
            "TP53": "tp53_scores.csv",
            "CHEK2": "chek2_scores.csv",
            "DNMT3A": "dnmt3a_garcia_2025_wt.csv",
        }
        dms_path = dms_dir / DMS_FILES[gene_symbol]

        try:
            with open(str(dms_path), "r") as fh:
                dms_rows = list(csv.DictReader(fh))
            gene_obj = get_gene(gene_symbol)
            dms_mapped = mapper_func(dms_rows, gene_obj)
        except Exception as e:
            log.warning(f"  Failed to load DMS for {gene_symbol}: {e}")
            continue

        if not dms_mapped:
            log.warning(f"  No DMS mappings for {gene_symbol}")
            continue

        # Build DMS lookup: variant_key -> DMS score
        # dms_mapped is a list of GenomicVariant dataclass objects
        dms_lookup = {}
        for gv in dms_mapped:
            key = f"{gv.chrom}:{gv.pos}:{gv.ref}>{gv.alt}"
            if gv.dms_score is not None:
                dms_lookup[key] = float(gv.dms_score)

        log.info(f"  DMS variants: {len(dms_lookup)}")

        # Get Evo2 scores for DMS-matched variants
        matched_keys = [k for k in dms_lookup if k in evo2_deltas]
        log.info(f"  Evo2-DMS matched: {len(matched_keys)}")

        if len(matched_keys) < 20:
            log.warning(f"  Too few matched variants")
            continue

        # Evo2 vs DMS
        dms_arr = np.array([dms_lookup[k] for k in matched_keys])
        evo2_arr = np.array([evo2_deltas[k] for k in matched_keys])
        rho_evo2, p_evo2 = scipy_stats.spearmanr(evo2_arr, dms_arr)
        log.info(f"  Evo2 vs DMS: rho={rho_evo2:.4f}, p={p_evo2:.2e}")

        gene_result = {
            "gene": gene_symbol,
            "n_dms": len(dms_lookup),
            "n_matched": len(matched_keys),
            "evo2_rho": float(rho_evo2),
            "evo2_p": float(p_evo2),
            "tool_rhos": {},
        }

        # Load other tool scores for same variants
        variant_key_set = set(matched_keys)
        tool_scores = load_tool_scores_for_gene(gene_symbol, variant_key_set)

        for tool_name, scores_dict in tool_scores.items():
            shared = [k for k in matched_keys if k in scores_dict]
            if len(shared) < 20:
                continue
            tool_arr = np.array([scores_dict[k] for k in shared])
            dms_sub = np.array([dms_lookup[k] for k in shared])
            rho, p = scipy_stats.spearmanr(tool_arr, dms_sub)
            gene_result["tool_rhos"][tool_name] = {
                "rho": float(rho),
                "p": float(p),
                "n": len(shared),
            }
            log.info(f"  {tool_name} vs DMS: rho={rho:.4f}, p={p:.2e} (n={len(shared)})")

        all_results.append(gene_result)

    # Save
    out_dir = OUT_DIR / "dms_multitool"
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "dms_multitool_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    log.info(f"\nSaved to {out_dir}")
    return all_results


# =============================================================================
# Analysis 5: Transversion control
# =============================================================================

def run_transversion_control():
    """Compare radiation-oxidative vs other transversions within each gene."""
    log.info("=" * 80)
    log.info("ANALYSIS 5: Transversion control")
    log.info("=" * 80)

    all_genes = [g.symbol for g in GENES.values()]
    results = []

    for gene_symbol in all_genes:
        log.info(f"\n--- {gene_symbol} ---")
        evo2_deltas = load_evo2_deltas(gene_symbol)
        if not evo2_deltas:
            continue

        # Classify variants by mutation type
        rad_oxidative = []  # C>A, G>T
        other_transversion = []  # A>C, A>T, G>C, T>A, T>G, C>G (excluding C>A, G>T)
        transitions = []  # C>T, T>C, A>G, G>A

        RAD_PAIRS = {("C", "A"), ("G", "T")}
        TRANSITION_PAIRS = {("C", "T"), ("T", "C"), ("A", "G"), ("G", "A")}
        OTHER_TV_PAIRS = {("A", "C"), ("A", "T"), ("G", "C"), ("T", "A"), ("T", "G"), ("C", "G")}

        for key, delta in evo2_deltas.items():
            parts = key.split(":")
            ref_alt = parts[2].split(">")
            if len(ref_alt) != 2 or len(ref_alt[0]) != 1 or len(ref_alt[1]) != 1:
                continue
            ref, alt = ref_alt
            pair = (ref, alt)
            if pair in RAD_PAIRS:
                rad_oxidative.append(delta)
            elif pair in OTHER_TV_PAIRS:
                other_transversion.append(delta)
            elif pair in TRANSITION_PAIRS:
                transitions.append(delta)

        rad_arr = np.array(rad_oxidative)
        otv_arr = np.array(other_transversion)
        trans_arr = np.array(transitions)

        # Test: radiation-oxidative vs other transversions
        if len(rad_arr) > 10 and len(otv_arr) > 10:
            u_stat, p_val = scipy_stats.mannwhitneyu(rad_arr, otv_arr, alternative="less")
            # rank-biserial correlation as effect size
            n1, n2 = len(rad_arr), len(otv_arr)
            rbc = 1 - (2 * u_stat) / (n1 * n2)

            log.info(f"  Rad-oxidative (n={len(rad_arr)}): mean={np.mean(rad_arr):.6f}")
            log.info(f"  Other transversion (n={len(otv_arr)}): mean={np.mean(otv_arr):.6f}")
            log.info(f"  Transitions (n={len(trans_arr)}): mean={np.mean(trans_arr):.6f}")
            log.info(f"  Rad vs other_TV: U={u_stat:.0f}, p={p_val:.2e}, rbc={rbc:.4f}")

            results.append({
                "gene": gene_symbol,
                "rad_oxidative_n": len(rad_arr),
                "rad_oxidative_mean": float(np.mean(rad_arr)),
                "rad_oxidative_std": float(np.std(rad_arr)),
                "other_transversion_n": len(otv_arr),
                "other_transversion_mean": float(np.mean(otv_arr)),
                "other_transversion_std": float(np.std(otv_arr)),
                "transition_n": len(trans_arr),
                "transition_mean": float(np.mean(trans_arr)),
                "transition_std": float(np.std(trans_arr)),
                "rad_vs_otherTV_U": float(u_stat),
                "rad_vs_otherTV_p": float(p_val),
                "rad_vs_otherTV_p_bonferroni": float(min(p_val * len(all_genes), 1.0)),
                "rad_vs_otherTV_rank_biserial": float(rbc),
                "rad_vs_transition_p": float(scipy_stats.mannwhitneyu(
                    rad_arr, trans_arr, alternative="less")[1]),
            })

    # Save
    out_dir = OUT_DIR / "transversion_control"
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "transversion_control_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    pd.DataFrame(results).to_csv(out_dir / "transversion_control_summary.csv", index=False)
    log.info(f"\nSaved to {out_dir}")
    return results


# =============================================================================
# Analysis 6: Ensemble model
# =============================================================================

def run_ensemble():
    """Evo2 + CADD logistic regression ensemble, 5-fold stratified CV."""
    log.info("=" * 80)
    log.info("ANALYSIS 6: Ensemble model (Evo2 + CADD)")
    log.info("=" * 80)

    all_results = []

    for gene_symbol in BENCHMARK_GENES:
        log.info(f"\n--- {gene_symbol} ---")

        evo2_deltas = load_evo2_deltas(gene_symbol)
        if not evo2_deltas:
            continue

        plp_keys, blb_keys = get_clinvar_variants(gene_symbol)
        all_keys = plp_keys + blb_keys
        all_labels = [1] * len(plp_keys) + [0] * len(blb_keys)

        # Load CADD scores
        variant_key_set = set(all_keys)
        tool_scores = load_tool_scores_for_gene(gene_symbol, variant_key_set)

        if "CADD" not in tool_scores:
            log.warning(f"  No CADD scores for {gene_symbol}")
            continue

        # Get variants with both Evo2 and CADD
        matched_keys = []
        matched_labels = []
        matched_evo2 = []
        matched_cadd = []

        for key, label in zip(all_keys, all_labels):
            if key in evo2_deltas and key in tool_scores["CADD"]:
                matched_keys.append(key)
                matched_labels.append(label)
                matched_evo2.append(-evo2_deltas[key])  # negate
                matched_cadd.append(tool_scores["CADD"][key])

        if len(matched_keys) < 20:
            continue

        labels = np.array(matched_labels)
        X_evo2 = np.array(matched_evo2).reshape(-1, 1)
        X_cadd = np.array(matched_cadd).reshape(-1, 1)
        X_both = np.column_stack([matched_evo2, matched_cadd])

        # 5-fold stratified CV
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        aurocs = {"evo2_only": [], "cadd_only": [], "ensemble": []}

        for train_idx, test_idx in skf.split(X_both, labels):
            # Evo2 only
            aurocs["evo2_only"].append(roc_auc_score(labels[test_idx], X_evo2[test_idx].ravel()))
            # CADD only
            aurocs["cadd_only"].append(roc_auc_score(labels[test_idx], X_cadd[test_idx].ravel()))
            # Ensemble
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_both[train_idx])
            X_test = scaler.transform(X_both[test_idx])
            lr = LogisticRegression(random_state=42, max_iter=1000)
            lr.fit(X_train, labels[train_idx])
            probs = lr.predict_proba(X_test)[:, 1]
            aurocs["ensemble"].append(roc_auc_score(labels[test_idx], probs))

        result = {
            "gene": gene_symbol,
            "n_variants": len(matched_keys),
            "n_plp": int(sum(labels)),
            "n_blb": int(sum(1 - labels)),
        }
        for model_name, vals in aurocs.items():
            result[f"{model_name}_mean_auroc"] = float(np.mean(vals))
            result[f"{model_name}_std_auroc"] = float(np.std(vals))
            result[f"{model_name}_folds"] = [float(v) for v in vals]
            log.info(f"  {model_name}: {np.mean(vals):.4f} ± {np.std(vals):.4f}")

        # Improvement
        evo2_mean = np.mean(aurocs["evo2_only"])
        cadd_mean = np.mean(aurocs["cadd_only"])
        ensemble_mean = np.mean(aurocs["ensemble"])
        best_single = max(evo2_mean, cadd_mean)
        result["improvement_over_best_single"] = float(ensemble_mean - best_single)
        result["improvement_pct"] = float((ensemble_mean - best_single) / best_single * 100)
        log.info(f"  Ensemble improvement: {result['improvement_pct']:.2f}% over best single tool")

        all_results.append(result)

    # Save
    out_dir = OUT_DIR / "ensemble"
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "ensemble_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    log.info(f"\nSaved to {out_dir}")
    return all_results


# =============================================================================
# Analysis 7: MPRA multi-tool comparison
# =============================================================================

def run_mpra_multitool():
    """Compute Spearman rho vs TERT MPRA for all tools."""
    log.info("=" * 80)
    log.info("ANALYSIS 7: TERT MPRA multi-tool comparison")
    log.info("=" * 80)

    gene_symbol = "TERT"
    gene = get_gene(gene_symbol)

    # Load Evo2 scores
    evo2_deltas = load_evo2_deltas(gene_symbol)

    # Load MPRA data — mapper needs (mpra_rows, reference_start, gene)
    TERT_MPRA_REFERENCE_START = 1294988
    mpra_path = MPRA_DIR / "tert" / "tert_gbm_scores.csv"
    try:
        with open(str(mpra_path), "r") as fh:
            mpra_rows = list(csv.DictReader(fh))
        mpra_mapped = map_tert_mpra(mpra_rows, TERT_MPRA_REFERENCE_START, gene)
    except Exception as e:
        log.error(f"  Failed to load MPRA: {e}")
        return None
    if not mpra_mapped:
        log.error("  No MPRA mappings")
        return None

    # Build MPRA lookup — GenomicVariant dataclass objects
    mpra_lookup = {}
    for gv in mpra_mapped:
        key = f"{gv.chrom}:{gv.pos}:{gv.ref}>{gv.alt}"
        if gv.dms_score is not None:
            mpra_lookup[key] = float(gv.dms_score)

    log.info(f"  MPRA variants: {len(mpra_lookup)}")

    # Get matched variants
    matched_keys = [k for k in mpra_lookup if k in evo2_deltas]
    log.info(f"  Evo2-MPRA matched: {len(matched_keys)}")

    if len(matched_keys) < 20:
        log.error("  Too few matched")
        return None

    mpra_arr = np.array([mpra_lookup[k] for k in matched_keys])
    evo2_arr = np.array([evo2_deltas[k] for k in matched_keys])

    rho_evo2, p_evo2 = scipy_stats.spearmanr(evo2_arr, mpra_arr)
    log.info(f"  Evo2 vs MPRA: rho={rho_evo2:.4f}, p={p_evo2:.2e}, n={len(matched_keys)}")

    result = {
        "gene": gene_symbol,
        "n_mpra": len(mpra_lookup),
        "n_matched": len(matched_keys),
        "evo2_rho": float(rho_evo2),
        "evo2_p": float(p_evo2),
        "tool_rhos": {},
    }

    # Load other tool scores
    variant_key_set = set(matched_keys)
    tool_scores = load_tool_scores_for_gene(gene_symbol, variant_key_set)

    for tool_name, scores_dict in tool_scores.items():
        shared = [k for k in matched_keys if k in scores_dict]
        if len(shared) < 20:
            result["tool_rhos"][tool_name] = {"rho": None, "n": len(shared), "note": "too few"}
            continue
        tool_arr = np.array([scores_dict[k] for k in shared])
        mpra_sub = np.array([mpra_lookup[k] for k in shared])
        rho, p = scipy_stats.spearmanr(tool_arr, mpra_sub)
        result["tool_rhos"][tool_name] = {
            "rho": float(rho),
            "p": float(p),
            "n": len(shared),
        }
        log.info(f"  {tool_name} vs MPRA: rho={rho:.4f}, p={p:.2e} (n={len(shared)})")

    # Save
    out_dir = OUT_DIR / "mpra_multitool"
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "mpra_multitool_results.json", "w") as f:
        json.dump(result, f, indent=2, default=str)
    log.info(f"\nSaved to {out_dir}")
    return result


# =============================================================================
# Analysis 8: Regenerate missing calibration files
# =============================================================================

def run_regenerate_calibration():
    """Regenerate DMS calibration JSONs for BRCA1 and CHEK2."""
    log.info("=" * 80)
    log.info("ANALYSIS 8: Regenerate calibration files")
    log.info("=" * 80)

    # This reuses DMS multi-tool analysis but saves in calibration format
    DMS_LOADERS = {
        "BRCA1": ("brca1", map_brca1_dms, "brca1_findlay_sge.csv"),
        "TP53": ("tp53", map_tp53_dms, "tp53_scores.csv"),
        "CHEK2": ("chek2", map_chek2_dms, "chek2_scores.csv"),
        "DNMT3A": ("dnmt3a", map_dnmt3a_garcia_dms, "dnmt3a_garcia_2025_wt.csv"),
    }

    for gene_symbol in DMS_GENES:
        log.info(f"\n--- {gene_symbol} ---")
        evo2_deltas = load_evo2_deltas(gene_symbol)
        if not evo2_deltas:
            continue

        subdir, mapper_func, filename = DMS_LOADERS[gene_symbol]
        dms_path = DMS_DIR / subdir / filename

        try:
            with open(str(dms_path), "r") as fh:
                dms_rows = list(csv.DictReader(fh))
            dms_mapped = mapper_func(dms_rows, get_gene(gene_symbol))
        except Exception as e:
            log.warning(f"  DMS load failed: {e}")
            continue

        # Build lookup — GenomicVariant dataclass objects
        dms_lookup = {}
        for gv in dms_mapped:
            key = f"{gv.chrom}:{gv.pos}:{gv.ref}>{gv.alt}"
            if gv.dms_score is not None:
                dms_lookup[key] = float(gv.dms_score)

        matched_keys = [k for k in dms_lookup if k in evo2_deltas]
        if len(matched_keys) < 10:
            continue

        dms_arr = np.array([dms_lookup[k] for k in matched_keys])
        evo2_arr = np.array([evo2_deltas[k] for k in matched_keys])
        rho, p = scipy_stats.spearmanr(evo2_arr, dms_arr)

        cal_result = {
            "gene": gene_symbol,
            "n_dms_total": len(dms_lookup),
            "n_evo2_matched": len(matched_keys),
            "spearman_rho": float(rho),
            "spearman_p": float(p),
            "evo2_mean_delta": float(np.mean(evo2_arr)),
            "evo2_std_delta": float(np.std(evo2_arr)),
            "dms_mean": float(np.mean(dms_arr)),
            "dms_std": float(np.std(dms_arr)),
        }

        out_dir = RESULTS_DIR / "calibration"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{gene_symbol}_calibration.json"
        with open(out_path, "w") as f:
            json.dump(cal_result, f, indent=2)
        log.info(f"  Saved: {out_path} (rho={rho:.4f}, p={p:.2e}, n={len(matched_keys)})")

    return True


# =============================================================================
# Analysis 9: Regenerate indel AUROCs
# =============================================================================

def run_regenerate_indels():
    """Compute frameshift vs B/LB ClinVar AUROCs for all genes."""
    log.info("=" * 80)
    log.info("ANALYSIS 9: Regenerate indel AUROCs")
    log.info("=" * 80)

    results = []

    for gene_symbol in [g.symbol for g in GENES.values()]:
        log.info(f"\n--- {gene_symbol} ---")

        # Load ALL Evo2 results (not just SNVs)
        ckpt_dir = RESULTS_DIR / gene_symbol / "w8192"
        checkpoint = ScoringCheckpoint(str(ckpt_dir))
        all_results_raw = checkpoint.load_results()

        if not all_results_raw:
            continue

        # Separate SNVs and indels
        snv_deltas = {}
        indel_deltas = {}
        for r in all_results_raw:
            key = f"{r['chrom']}:{r['pos']}:{r['ref']}>{r['alt']}"
            if len(r["ref"]) == 1 and len(r["alt"]) == 1:
                snv_deltas[key] = r["delta"]
            else:
                indel_deltas[key] = r["delta"]

        # Get ClinVar indels
        entries_plp = parse_clinvar_vcf(
            str(CLINVAR_VCF), genes={gene_symbol},
            min_stars=1, significance_filter={"P/LP"},
        )
        entries_blb = parse_clinvar_vcf(
            str(CLINVAR_VCF), genes={gene_symbol},
            min_stars=1, significance_filter={"B/LB"},
        )

        plp_indel_keys = []
        blb_indel_keys = []
        plp_fs_keys = []  # frameshift
        plp_if_keys = []  # inframe

        for e in clinvar_to_variants(entries_plp):
            if not e.is_snv and e.key() in indel_deltas:
                plp_indel_keys.append(e.key())
                # Classify as frameshift vs inframe
                ref_len = len(e.ref)
                alt_len = len(e.alt)
                if (alt_len - ref_len) % 3 != 0:
                    plp_fs_keys.append(e.key())
                else:
                    plp_if_keys.append(e.key())

        for e in clinvar_to_variants(entries_blb):
            if not e.is_snv and e.key() in indel_deltas:
                blb_indel_keys.append(e.key())

        # Also use SNV B/LB as negative controls if indel B/LB too few
        blb_snv_keys = []
        for e in clinvar_to_variants(entries_blb):
            if e.is_snv and e.key() in snv_deltas:
                blb_snv_keys.append(e.key())

        log.info(f"  Indels: {len(plp_indel_keys)} P/LP ({len(plp_fs_keys)} fs, {len(plp_if_keys)} if), {len(blb_indel_keys)} B/LB")

        gene_result = {
            "gene": gene_symbol,
            "n_indel_plp": len(plp_indel_keys),
            "n_indel_blb": len(blb_indel_keys),
            "n_frameshift_plp": len(plp_fs_keys),
            "n_inframe_plp": len(plp_if_keys),
        }

        # Frameshift AUROC (if enough frameshifts and B/LB)
        if len(plp_fs_keys) >= 3 and len(blb_indel_keys) >= 3:
            labels = np.array([1] * len(plp_fs_keys) + [0] * len(blb_indel_keys))
            scores = np.array(
                [-indel_deltas[k] for k in plp_fs_keys] +
                [-indel_deltas[k] for k in blb_indel_keys]
            )
            try:
                auroc = roc_auc_score(labels, scores)
                gene_result["frameshift_auroc"] = float(auroc)
                gene_result["frameshift_n"] = len(labels)
                log.info(f"  Frameshift AUROC: {auroc:.4f} (n={len(labels)})")
            except ValueError:
                pass

        # Mean deltas
        if plp_fs_keys:
            gene_result["frameshift_mean_delta"] = float(np.mean([indel_deltas[k] for k in plp_fs_keys]))
        if plp_if_keys:
            gene_result["inframe_mean_delta"] = float(np.mean([indel_deltas[k] for k in plp_if_keys]))

        # Total indel stats
        gene_result["total_indels_scored"] = len(indel_deltas)
        gene_result["total_snvs_scored"] = len(snv_deltas)

        results.append(gene_result)

    # Save
    out_dir = OUT_DIR / "indels"
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "indel_aurocs_all_genes.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    pd.DataFrame(results).to_csv(out_dir / "indel_summary.csv", index=False)
    log.info(f"\nSaved to {out_dir}")
    return results


# =============================================================================
# Analysis 10: Exact variant count
# =============================================================================

def run_variant_count():
    """Compute exact total variant count across all 10 genes."""
    log.info("=" * 80)
    log.info("ANALYSIS 10: Exact variant count")
    log.info("=" * 80)

    total = 0
    total_snvs = 0
    total_indels = 0
    gene_counts = []

    for gene_symbol in [g.symbol for g in GENES.values()]:
        ckpt_dir = RESULTS_DIR / gene_symbol / "w8192"
        checkpoint = ScoringCheckpoint(str(ckpt_dir))
        results = checkpoint.load_results()

        n_snv = 0
        n_indel = 0
        for r in results:
            if len(r["ref"]) == 1 and len(r["alt"]) == 1:
                n_snv += 1
            else:
                n_indel += 1

        n_total = len(results)
        total += n_total
        total_snvs += n_snv
        total_indels += n_indel

        gene_counts.append({
            "gene": gene_symbol,
            "total": n_total,
            "snvs": n_snv,
            "indels": n_indel,
        })
        log.info(f"  {gene_symbol}: {n_total} total ({n_snv} SNVs, {n_indel} indels)")

    log.info(f"\n  TOTAL: {total} variants ({total_snvs} SNVs, {total_indels} indels)")

    result = {
        "total_variants": total,
        "total_snvs": total_snvs,
        "total_indels": total_indels,
        "per_gene": gene_counts,
    }

    out_dir = OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "variant_counts.json", "w") as f:
        json.dump(result, f, indent=2)
    log.info(f"\nSaved to {out_dir / 'variant_counts.json'}")
    return result


# =============================================================================
# Main
# =============================================================================

ANALYSES = {
    "matched_benchmark": run_matched_benchmark,
    "dms_multitool": run_dms_multitool,
    "transversion_control": run_transversion_control,
    "ensemble": run_ensemble,
    "mpra_multitool": run_mpra_multitool,
    "regenerate_calibration": run_regenerate_calibration,
    "regenerate_indels": run_regenerate_indels,
    "variant_count": run_variant_count,
}


def main():
    parser = argparse.ArgumentParser(description="Pre-draft manuscript analyses")
    parser.add_argument(
        "--analysis",
        type=str,
        default="all",
        choices=list(ANALYSES.keys()) + ["all"],
        help="Which analysis to run",
    )
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.analysis == "all":
        for name, func in ANALYSES.items():
            try:
                func()
            except Exception as e:
                log.error(f"Analysis '{name}' failed: {e}")
                import traceback
                traceback.print_exc()
    else:
        ANALYSES[args.analysis]()


if __name__ == "__main__":
    main()

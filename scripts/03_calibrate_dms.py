#!/usr/bin/env python3
"""
DMS calibration for control genes (Phase 2).

For each control gene (BRCA1, TP53, CHEK2, DNMT3A):
1. Load DMS fitness scores from MaveDB/supplements
2. Load Evo2 scoring results
3. Compute Spearman rho between Evo2 delta and DMS fitness
4. Compute gene-specific constraint landscape (all possible SNVs)
5. Calibrate gene-specific thresholds via Pejaver 2022 LR

Go/no-go gate: BRCA1 AUROC ≥ 0.85

Produces: Figure 3, Figure 5 (control panels), Supp Fig S6

Usage:
  python 03_calibrate_dms.py --gene BRCA1
  python 03_calibrate_dms.py --gene all
"""

import argparse
import csv
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from utils.config import DMS_DIR, MPRA_DIR, RESULTS_DIR
from utils.gene_coordinates import CONTROL_GENES, GENES, get_gene
from utils.hgvs_mapping import (
    build_cds_map,
    map_brca1_dms,
    map_chek2_dms,
    map_dnmt3a_dms,
    map_dnmt3a_garcia_dms,
    map_tp53_dms,
    map_tert_mpra,
)
from utils.sequence_utils import ScoringCheckpoint
from utils.benchmarking import (
    calibrate_gene_thresholds,
    compute_auroc,
    compute_auprc,
    compute_score_landscape,
    run_sanity_checks,
)

try:
    from scipy import stats as scipy_stats
except ImportError:
    scipy_stats = None

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)


# =============================================================================
# DMS Data Loaders
# =============================================================================

def load_dms_mapped(gene_symbol: str, dms_dir: Path) -> list:
    """
    Load DMS data and map to genomic coordinates via HGVS mapping.

    Returns list of GenomicVariant objects with dms_score set.
    """
    gene = get_gene(gene_symbol)

    if gene_symbol == "BRCA1":
        filepath = dms_dir / "brca1" / "brca1_findlay_sge.csv"
        if not filepath.exists():
            log.warning(f"BRCA1 DMS file not found: {filepath}")
            return []
        with open(filepath) as f:
            rows = list(csv.DictReader(f))
        log.info(f"  BRCA1 DMS: {len(rows)} rows loaded")
        variants = map_brca1_dms(rows, gene)
        log.info(f"  BRCA1 DMS: {len(variants)} SNVs mapped to genomic coords")
        return variants

    elif gene_symbol == "TP53":
        filepath = dms_dir / "tp53" / "tp53_scores.csv"
        if not filepath.exists():
            log.warning(f"TP53 DMS file not found: {filepath}")
            return []
        with open(filepath) as f:
            rows = list(csv.DictReader(f))
        log.info(f"  TP53 DMS: {len(rows)} rows loaded")
        # TP53 has protein-level HGVS only — returns multiple candidates per AA change
        variants = map_tp53_dms(rows, gene)
        log.info(f"  TP53 DMS: {len(variants)} candidate SNVs mapped")
        return variants

    elif gene_symbol == "CHEK2":
        filepath = dms_dir / "chek2" / "chek2_scores.csv"
        if not filepath.exists():
            log.warning(f"CHEK2 DMS file not found: {filepath}")
            return []
        with open(filepath) as f:
            rows = list(csv.DictReader(f))
        log.info(f"  CHEK2 DMS: {len(rows)} rows loaded")
        variants = map_chek2_dms(rows, gene, score_column="RCS_this_SNV")
        log.info(f"  CHEK2 DMS: {len(variants)} SNVs mapped to genomic coords")
        return variants

    elif gene_symbol == "DNMT3A":
        # Prefer Garcia 2025 methylation activity (2,036 variants) over Garcia 2025 (2,036) is primary; Huang 2022 stability (254) is fallback
        garcia_path = dms_dir / "dnmt3a" / "dnmt3a_garcia_2025_wt.csv"
        huang_path = dms_dir / "dnmt3a" / "dnmt3a_scores.csv"
        if garcia_path.exists():
            with open(garcia_path) as f:
                rows = list(csv.DictReader(f))
            log.info(f"  DNMT3A DMS: {len(rows)} rows loaded (Garcia 2025 methylation activity)")
            variants = map_dnmt3a_garcia_dms(rows, gene, score_column="score")
            log.info(f"  DNMT3A DMS: {len(variants)} candidate SNVs mapped to genomic coords")
            return variants
        elif huang_path.exists():
            with open(huang_path) as f:
                rows = list(csv.DictReader(f))
            log.info(f"  DNMT3A DMS: {len(rows)} rows loaded (Huang 2022 stability, fallback — Garcia 2025 preferred)")
            variants = map_dnmt3a_dms(rows, gene)
            log.info(f"  DNMT3A DMS: {len(variants)} SNVs mapped to genomic coords")
            return variants
        else:
            log.warning(f"DNMT3A DMS file not found: tried {garcia_path} and {huang_path}")
            return []

    else:
        log.warning(f"  No DMS loader for {gene_symbol}")
        return []


# =============================================================================
# Calibration
# =============================================================================

def calibrate_gene(gene_symbol: str, window_size: int = 8192):
    """Run full DMS calibration for one control gene."""
    gene = get_gene(gene_symbol)
    log.info(f"\n{'='*60}")
    log.info(f"Calibrating {gene_symbol}")
    log.info(f"  DMS source: {gene.dms_source[:80]}...")
    log.info(f"{'='*60}")

    # Load DMS data with HGVS-to-genomic mapping
    dms_variants = load_dms_mapped(gene_symbol, DMS_DIR)
    if not dms_variants:
        log.warning(f"No mapped DMS data for {gene_symbol}")

    # Load Evo2 scoring results
    ckpt_dir = RESULTS_DIR / gene_symbol / f"w{window_size}"
    checkpoint = ScoringCheckpoint(str(ckpt_dir))
    evo2_results = checkpoint.load_results()
    log.info(f"  Evo2 results: {len(evo2_results)} variants scored")

    if not evo2_results:
        log.error(f"No Evo2 results for {gene_symbol}")
        return None

    # Build Evo2 delta lookup by variant key
    evo2_deltas = {}
    for r in evo2_results:
        key = f"{r['chrom']}:{r['pos']}:{r['ref']}>{r['alt']}"
        evo2_deltas[key] = r["delta"]

    # Match DMS scores to Evo2 deltas via genomic coordinates
    dms_matched = []
    if dms_variants:
        for dv in dms_variants:
            key = f"{dv.chrom}:{dv.pos}:{dv.ref}>{dv.alt}"
            if key in evo2_deltas:
                dms_matched.append({
                    "key": key,
                    "dms_score": dv.dms_score,
                    "evo2_delta": evo2_deltas[key],
                    "hgvs": dv.hgvs,
                })
        log.info(f"  DMS-Evo2 matches: {len(dms_matched)} / {len(dms_variants)} DMS variants")

    # Compute Spearman rho between DMS and Evo2
    spearman_rho = None
    spearman_p = None
    if len(dms_matched) >= 10 and scipy_stats is not None:
        dms_scores = np.array([m["dms_score"] for m in dms_matched])
        evo2_scores = np.array([m["evo2_delta"] for m in dms_matched])
        spearman_rho, spearman_p = scipy_stats.spearmanr(dms_scores, evo2_scores)
        log.info(f"  DMS Spearman rho: {spearman_rho:.4f} (p={spearman_p:.2e}, n={len(dms_matched)})")

    # ClinVar P/LP vs B/LB validation
    plp_deltas = [r["delta"] for r in evo2_results if r.get("clinvar_class") == "P/LP"]
    blb_deltas = [r["delta"] for r in evo2_results if r.get("clinvar_class") == "B/LB"]

    auroc = None
    thresholds = {}
    if plp_deltas and blb_deltas:
        labels = np.array([1] * len(plp_deltas) + [0] * len(blb_deltas))
        scores = np.array(plp_deltas + blb_deltas)
        auroc, _, _, _ = compute_auroc(labels, -scores)
        log.info(f"  ClinVar AUROC: {auroc:.4f} ({len(plp_deltas)} P/LP, {len(blb_deltas)} B/LB)")

        # Calibrate thresholds
        thresholds = calibrate_gene_thresholds(
            pathogenic_scores=-np.array(plp_deltas),
            benign_scores=-np.array(blb_deltas),
        )
        log.info(f"  Calibrated thresholds: {thresholds}")
    else:
        log.warning(f"  No ClinVar P/LP or B/LB variants to validate against")

    # Score landscape (all possible SNVs)
    all_snv_deltas = [r["delta"] for r in evo2_results if r.get("region_type") == "coding"]
    if all_snv_deltas:
        landscape = compute_score_landscape(
            all_snv_scores=-np.array(all_snv_deltas),
            pathogenic_scores=-np.array(plp_deltas) if plp_deltas else None,
            benign_scores=-np.array(blb_deltas) if blb_deltas else None,
        )
    else:
        landscape = {}

    # Sanity checks
    checks = run_sanity_checks(evo2_results, gene_symbol)
    log.info(f"  Sanity checks: {checks}")

    # Save calibration results
    result = {
        "gene": gene_symbol,
        "window_size": window_size,
        "n_evo2_scored": len(evo2_results),
        "n_plp": len(plp_deltas),
        "n_blb": len(blb_deltas),
        "n_dms_mapped": len(dms_variants),
        "n_dms_matched": len(dms_matched),
        "spearman_rho": spearman_rho,
        "spearman_p": spearman_p,
        "auroc": auroc,
        "thresholds": thresholds,
        "sanity_checks": checks,
    }

    out_dir = RESULTS_DIR / gene_symbol / "calibration"
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "calibration_result.json", 'w') as f:
        json.dump(result, f, indent=2, default=str)

    # Save DMS-Evo2 matched pairs for downstream plotting
    if dms_matched:
        with open(out_dir / "dms_evo2_matched.json", 'w') as f:
            json.dump(dms_matched, f, indent=2)
        log.info(f"  Saved {len(dms_matched)} matched pairs to {out_dir / 'dms_evo2_matched.json'}")

    if landscape:
        with open(out_dir / "score_landscape.json", 'w') as f:
            json.dump(landscape, f, indent=2)

    return result


def main():
    parser = argparse.ArgumentParser(description="DMS calibration for control genes")
    parser.add_argument("--gene", type=str, default="all",
                        help="Gene symbol or 'all' (default: all 4 controls)")
    parser.add_argument("--window-size", type=int, default=8192)
    args = parser.parse_args()

    control_syms = [g.symbol for g in CONTROL_GENES]
    genes = control_syms if args.gene == "all" else [args.gene.upper()]

    log.info(f"DMS Calibration — genes: {genes}")

    results = []
    for gene_sym in genes:
        if gene_sym not in control_syms:
            log.warning(f"{gene_sym} is not a control gene, skipping")
            continue
        result = calibrate_gene(gene_sym, args.window_size)
        if result:
            results.append(result)

    # Go/no-go gate
    log.info(f"\n{'='*60}")
    log.info("GO/NO-GO GATE")
    log.info(f"{'='*60}")
    brca1_result = next((r for r in results if r["gene"] == "BRCA1"), None)
    if brca1_result and brca1_result["auroc"] is not None:
        if brca1_result["auroc"] >= 0.85:
            log.info(f"  BRCA1 AUROC = {brca1_result['auroc']:.4f} >= 0.85 → GO")
        else:
            log.warning(
                f"  BRCA1 AUROC = {brca1_result['auroc']:.4f} < 0.85 → NO-GO "
                f"(debug before proceeding)"
            )
    else:
        log.warning("  BRCA1 calibration not completed — cannot evaluate gate")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
ClinVar independent validation (Phase 3).

Validates Evo2 VEP scores against ClinVar P/LP vs B/LB variants.
Separate from DMS calibration (Phase 2) to avoid circularity.

For each gene:
1. Load Evo2 scoring results (from Phase 2 / combined scoring)
2. Stratify by ClinVar review status (1-star, ≥2-star, ≥3-star)
3. Compute AUROC, AUPRC, precision-recall curves, calibration curves
4. Compare against benchmark tools (REVEL, CADD, AlphaMissense)
5. Compute orthogonality (pairwise Spearman) between Evo2 and other tools

Produces: Figure 4, Supp Fig S1, S4, S5

Usage:
  python 04_validate_clinvar.py --gene BRCA1
  python 04_validate_clinvar.py --gene all
  python 04_validate_clinvar.py --gene all --star-sensitivity
"""

import argparse
import csv
import gzip
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from utils.config import (
    ALPHAMISSENSE_TSV,
    CLINVAR_VCF,
    RESULTS_DIR,
    SHARED_DATA_DIR,
)
from utils.gene_coordinates import GENES, get_gene
from utils.clinvar_parser import parse_clinvar_vcf, clinvar_to_variants
from utils.sequence_utils import ScoringCheckpoint
from utils.benchmarking import (
    ToolComparison,
    compute_auroc,
    compute_auprc,
    compute_calibration,
    compute_orthogonality_matrix,
    interpret_orthogonality,
    precision_at_recall,
    recall_at_precision,
    save_comparison,
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

REVEL_PATH = SHARED_DATA_DIR / "benchmarks" / "revel" / "revel_genes.bed.gz"
CADD_PATH = SHARED_DATA_DIR / "benchmarks" / "cadd" / "whole_genome_SNVs.tsv.gz"
CADD_TBI = SHARED_DATA_DIR / "benchmarks" / "cadd" / "whole_genome_SNVs.tsv.gz.tbi"


# =============================================================================
# Benchmark Score Loaders
# =============================================================================

def load_alphamissense_scores(
    gene_chrom: str,
    gene_start: int,
    gene_end: int,
) -> Dict[str, float]:
    """
    Load AlphaMissense scores for a genomic region.

    Returns dict: variant_key -> AM pathogenicity score (0-1).
    """
    scores = {}
    path = str(ALPHAMISSENSE_TSV)

    try:
        import pysam
        tbx = pysam.TabixFile(path)
        # AM uses chr-prefixed coordinates, 1-based in the TSV
        chrom_query = gene_chrom if gene_chrom.startswith("chr") else "chr" + gene_chrom
        for row in tbx.fetch(chrom_query, gene_start, gene_end):
            fields = row.split("\t")
            # Columns: CHROM, POS, REF, ALT, genome, uniprot_id, transcript_id,
            #          protein_variant, am_pathogenicity, am_class
            chrom = fields[0] if fields[0].startswith("chr") else "chr" + fields[0]
            pos = int(fields[1]) - 1  # Convert to 0-based
            ref = fields[2]
            alt = fields[3]
            am_score = float(fields[8])
            key = "{}:{}:{}>{}".format(chrom, pos, ref, alt)
            scores[key] = am_score
        tbx.close()
    except Exception as e:
        log.warning("Could not load AlphaMissense: {}".format(e))

    return scores


def load_cadd_scores(
    gene_chrom: str,
    gene_start: int,
    gene_end: int,
) -> Dict[str, float]:
    """
    Load CADD v1.7 PHRED scores for a genomic region via tabix.

    Returns dict: variant_key -> CADD PHRED score.
    """
    scores = {}
    path = str(CADD_PATH)

    if not Path(path).exists():
        log.warning("CADD file not found: {}".format(path))
        return scores

    try:
        import pysam
        tbx = pysam.TabixFile(path)
        # CADD uses non-chr coordinates, 1-based
        chrom_query = gene_chrom.replace("chr", "")
        for row in tbx.fetch(chrom_query, gene_start, gene_end):
            fields = row.split("\t")
            # Columns: Chrom, Pos, Ref, Alt, RawScore, PHRED
            chrom = "chr" + fields[0]
            pos = int(fields[1]) - 1  # Convert to 0-based
            ref = fields[2]
            alt = fields[3]
            phred = float(fields[5])
            key = "{}:{}:{}>{}".format(chrom, pos, ref, alt)
            scores[key] = phred
        tbx.close()
    except Exception as e:
        log.warning("Could not load CADD: {}".format(e))

    return scores


def load_revel_scores(
    gene_chrom: str,
    gene_start: int,
    gene_end: int,
) -> Dict[str, float]:
    """
    Load REVEL scores for a genomic region via tabix.

    BED format: chrom, start(0-based), end, ref, alt, aa_ref, aa_alt, revel_score
    Returns dict: variant_key -> REVEL score (0-1).
    """
    scores = {}
    path = str(REVEL_PATH)

    if not Path(path).exists():
        log.warning("REVEL file not found: {}".format(path))
        return scores

    try:
        import pysam
        tbx = pysam.TabixFile(path)
        chrom_query = gene_chrom if gene_chrom.startswith("chr") else "chr" + gene_chrom
        for row in tbx.fetch(chrom_query, gene_start, gene_end):
            fields = row.split("\t")
            # BED: chrom, start(0-based), end, ref, alt, aa_ref, aa_alt, revel
            chrom = fields[0]
            pos = int(fields[1])  # already 0-based
            ref = fields[3]
            alt = fields[4]
            revel = float(fields[7])
            key = "{}:{}:{}>{}".format(chrom, pos, ref, alt)
            scores[key] = revel
        tbx.close()
    except Exception as e:
        log.warning("Could not load REVEL: {}".format(e))

    return scores


# =============================================================================
# Validation
# =============================================================================

def validate_gene(
    gene_symbol: str,
    window_size: int = 8192,
    star_sensitivity: bool = False,
) -> Optional[dict]:
    """Run full ClinVar validation for one gene."""
    gene = get_gene(gene_symbol)
    log.info("\n{}".format("=" * 60))
    log.info("ClinVar Validation: {}".format(gene_symbol))
    log.info("=" * 60)

    # Load Evo2 scoring results
    ckpt_dir = RESULTS_DIR / gene_symbol / "w{}".format(window_size)
    checkpoint = ScoringCheckpoint(str(ckpt_dir))
    evo2_results = checkpoint.load_results()
    log.info("  Evo2 results: {} variants scored".format(len(evo2_results)))

    if not evo2_results:
        log.warning("  No Evo2 results for {}".format(gene_symbol))
        return None

    # Build Evo2 delta lookup
    evo2_deltas = {}
    for r in evo2_results:
        key = "{}:{}:{}>{}".format(r["chrom"], r["pos"], r["ref"], r["alt"])
        evo2_deltas[key] = r["delta"]

    # Get ClinVar variants
    entries_plp = parse_clinvar_vcf(
        str(CLINVAR_VCF), genes={gene_symbol},
        min_stars=1, significance_filter={"P/LP"},
    )
    entries_blb = parse_clinvar_vcf(
        str(CLINVAR_VCF), genes={gene_symbol},
        min_stars=1, significance_filter={"B/LB"},
    )

    variants_plp = clinvar_to_variants(entries_plp)
    variants_blb = clinvar_to_variants(entries_blb)

    # Only SNVs
    variants_plp = [v for v in variants_plp if v.is_snv]
    variants_blb = [v for v in variants_blb if v.is_snv]

    log.info("  ClinVar SNVs: {} P/LP, {} B/LB (≥1 star)".format(
        len(variants_plp), len(variants_blb)))

    # Match to Evo2 results
    matched_plp = []
    matched_blb = []
    all_variant_keys = []

    for v in variants_plp:
        key = v.key()
        if key in evo2_deltas:
            matched_plp.append({
                "key": key, "delta": evo2_deltas[key],
                "stars": v.clinvar_stars,
            })
            all_variant_keys.append(key)

    for v in variants_blb:
        key = v.key()
        if key in evo2_deltas:
            matched_blb.append({
                "key": key, "delta": evo2_deltas[key],
                "stars": v.clinvar_stars,
            })
            all_variant_keys.append(key)

    log.info("  Matched to Evo2: {} P/LP, {} B/LB".format(
        len(matched_plp), len(matched_blb)))

    if len(matched_plp) < 5 or len(matched_blb) < 5:
        log.warning("  Too few matched variants for validation")
        return None

    # Star-level stratification
    star_results = {}
    for min_stars in [1, 2, 3]:
        plp_deltas = [m["delta"] for m in matched_plp if m["stars"] >= min_stars]
        blb_deltas = [m["delta"] for m in matched_blb if m["stars"] >= min_stars]

        if len(plp_deltas) < 5 or len(blb_deltas) < 5:
            continue

        labels = np.array([1] * len(plp_deltas) + [0] * len(blb_deltas))
        scores = -np.array(plp_deltas + blb_deltas)

        auroc, fpr, tpr, _ = compute_auroc(labels, scores)
        auprc, prec, recall, _ = compute_auprc(labels, scores)
        p_at_r90 = precision_at_recall(labels, scores, 0.9)
        r_at_p90 = recall_at_precision(labels, scores, 0.9)

        star_results["ge_{}_star".format(min_stars)] = {
            "auroc": float(auroc),
            "auprc": float(auprc),
            "n_plp": len(plp_deltas),
            "n_blb": len(blb_deltas),
            "precision_at_90_recall": float(p_at_r90),
            "recall_at_90_precision": float(r_at_p90),
            "mean_delta_plp": float(np.mean(plp_deltas)),
            "mean_delta_blb": float(np.mean(blb_deltas)),
        }
        log.info("  ≥{}★: AUROC={:.4f}, AUPRC={:.4f}, P/LP={}, B/LB={}".format(
            min_stars, auroc, auprc, len(plp_deltas), len(blb_deltas)))

    # Primary validation: ≥2 star
    primary_key = "ge_2_star"
    if primary_key not in star_results:
        primary_key = "ge_1_star"

    # Load benchmark tool scores
    log.info("  Loading benchmark scores...")
    tool_scores = {"Evo2": {}}

    # Build Evo2 score array for all matched variants
    for m in matched_plp + matched_blb:
        tool_scores["Evo2"][m["key"]] = -m["delta"]

    # AlphaMissense
    am_scores = load_alphamissense_scores(gene.chrom, gene.start, gene.end)
    if am_scores:
        tool_scores["AlphaMissense"] = am_scores
        n_am = sum(1 for k in all_variant_keys if k in am_scores)
        log.info("  AlphaMissense: {} scores loaded, {} matched".format(
            len(am_scores), n_am))

    # CADD
    cadd_scores = load_cadd_scores(gene.chrom, gene.start, gene.end)
    if cadd_scores:
        tool_scores["CADD"] = cadd_scores
        n_cadd = sum(1 for k in all_variant_keys if k in cadd_scores)
        log.info("  CADD: {} scores loaded, {} matched".format(
            len(cadd_scores), n_cadd))

    # REVEL (tabix-indexed BED, fast regional query)
    revel_scores = load_revel_scores(gene.chrom, gene.start, gene.end)
    if revel_scores:
        tool_scores["REVEL"] = revel_scores
        n_revel = sum(1 for k in all_variant_keys if k in revel_scores)
        log.info("  REVEL: {} scores loaded, {} matched".format(
            len(revel_scores), n_revel))

    # Compute tool comparison (AUROC for each tool on shared variants)
    all_labels = ([1] * len(matched_plp)) + ([0] * len(matched_blb))
    all_keys = [m["key"] for m in matched_plp] + [m["key"] for m in matched_blb]
    labels_arr = np.array(all_labels)

    tool_aurocs = {}
    tool_auprcs = {}
    for tool_name, scores_dict in tool_scores.items():
        # Build score array aligned with labels
        tool_vals = []
        tool_labels = []
        for key, label in zip(all_keys, all_labels):
            if key in scores_dict:
                tool_vals.append(scores_dict[key])
                tool_labels.append(label)

        if len(tool_vals) < 10 or len(set(tool_labels)) < 2:
            continue

        tool_labels_arr = np.array(tool_labels)
        tool_vals_arr = np.array(tool_vals)

        auroc, _, _, _ = compute_auroc(tool_labels_arr, tool_vals_arr)
        auprc, _, _, _ = compute_auprc(tool_labels_arr, tool_vals_arr)
        tool_aurocs[tool_name] = float(auroc)
        tool_auprcs[tool_name] = float(auprc)
        log.info("  {}: AUROC={:.4f}, AUPRC={:.4f} (n={})".format(
            tool_name, auroc, auprc, len(tool_vals)))

    # Orthogonality (Spearman correlations between tools)
    ortho_results = {}
    if len(tool_scores) >= 2 and scipy_stats is not None:
        # Build aligned score arrays for shared variants
        tool_names = sorted(tool_scores.keys())
        for i, t1 in enumerate(tool_names):
            for j, t2 in enumerate(tool_names):
                if j <= i:
                    continue
                shared = [k for k in all_keys
                          if k in tool_scores[t1] and k in tool_scores[t2]]
                if len(shared) < 10:
                    continue
                s1 = np.array([tool_scores[t1][k] for k in shared])
                s2 = np.array([tool_scores[t2][k] for k in shared])
                rho, p = scipy_stats.spearmanr(s1, s2)
                pair_key = "{}_vs_{}".format(t1, t2)
                ortho_results[pair_key] = {
                    "spearman_rho": float(rho),
                    "p_value": float(p),
                    "n_shared": len(shared),
                    "interpretation": interpret_orthogonality(rho),
                }
                log.info("  {} vs {}: rho={:.4f} ({})".format(
                    t1, t2, rho, interpret_orthogonality(rho)))

    # Calibration curve (Evo2 scores)
    plp_deltas_all = [m["delta"] for m in matched_plp]
    blb_deltas_all = [m["delta"] for m in matched_blb]
    labels_all = np.array([1] * len(plp_deltas_all) + [0] * len(blb_deltas_all))
    scores_all = -np.array(plp_deltas_all + blb_deltas_all)

    try:
        frac_pos, mean_pred = compute_calibration(labels_all, scores_all, n_bins=10)
        calibration = {
            "fraction_positive": frac_pos.tolist(),
            "mean_predicted": mean_pred.tolist(),
        }
    except Exception:
        calibration = {}

    # Assemble results
    result = {
        "gene": gene_symbol,
        "window_size": window_size,
        "n_evo2_scored": len(evo2_results),
        "n_plp_matched": len(matched_plp),
        "n_blb_matched": len(matched_blb),
        "star_stratification": star_results,
        "tool_aurocs": tool_aurocs,
        "tool_auprcs": tool_auprcs,
        "orthogonality": ortho_results,
        "calibration": calibration,
    }

    # Save
    out_dir = RESULTS_DIR / gene_symbol / "clinvar_validation"
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(out_dir / "validation_result.json", "w") as f:
        json.dump(result, f, indent=2, default=str)

    log.info("  Saved to {}".format(out_dir / "validation_result.json"))
    return result


def main():
    parser = argparse.ArgumentParser(description="ClinVar independent validation")
    parser.add_argument("--gene", type=str, default="all",
                        help="Gene symbol or 'all'")
    parser.add_argument("--window-size", type=int, default=8192,
                        help="Window size used for scoring")
    parser.add_argument("--star-sensitivity", action="store_true",
                        help="Run star-level sensitivity analysis")
    args = parser.parse_args()

    all_symbols = [g.symbol for g in GENES.values()]
    genes = all_symbols if args.gene == "all" else [args.gene.upper()]

    log.info("ClinVar Validation — genes: {}".format(genes))

    results = []
    for gene_sym in genes:
        result = validate_gene(
            gene_sym, args.window_size,
            star_sensitivity=args.star_sensitivity,
        )
        if result:
            results.append(result)

    # Summary
    log.info("\n{}".format("=" * 80))
    log.info("CLINVAR VALIDATION SUMMARY")
    log.info("=" * 80)
    log.info("{:<8} {:>8} {:>8} {:>8} {:>6} {:>6}".format(
        "Gene", "AUROC", "AUPRC", "P@R90", "P/LP", "B/LB"))
    log.info("-" * 50)

    for r in results:
        primary = r["star_stratification"].get("ge_2_star",
                  r["star_stratification"].get("ge_1_star", {}))
        if primary:
            log.info("{:<8} {:>8.4f} {:>8.4f} {:>8.4f} {:>6} {:>6}".format(
                r["gene"],
                primary["auroc"],
                primary["auprc"],
                primary["precision_at_90_recall"],
                primary["n_plp"],
                primary["n_blb"],
            ))

    # Tool comparison summary
    log.info("\nTool AUROC Comparison:")
    all_tools = set()
    for r in results:
        all_tools.update(r["tool_aurocs"].keys())

    if all_tools:
        header = "{:<8}".format("Gene")
        for t in sorted(all_tools):
            header += " {:>12}".format(t[:12])
        log.info(header)
        log.info("-" * (8 + 13 * len(all_tools)))

        for r in results:
            line = "{:<8}".format(r["gene"])
            for t in sorted(all_tools):
                auroc = r["tool_aurocs"].get(t)
                if auroc is not None:
                    line += " {:>12.4f}".format(auroc)
                else:
                    line += " {:>12}".format("N/A")
            log.info(line)

    # Save aggregate
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_DIR / "clinvar_validation_summary.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    log.info("\nAggregate saved to {}".format(
        RESULTS_DIR / "clinvar_validation_summary.json"))


if __name__ == "__main__":
    main()

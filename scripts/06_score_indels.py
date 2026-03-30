#!/usr/bin/env python3
"""
Indel scoring with Evo2 (Phase 4 supplement).

Scores in-frame and frameshift indels separately:
  - Frameshifts: trivially LOF; reported as sanity check (should score strongly negative)
  - In-frame indels: novel — Evo2 assesses impact on local sequence grammar

Produces: Supp Fig S3 data

Usage:
  python 06_score_indels.py --gene BRCA1 --window-size 8192
  python 06_score_indels.py --gene all
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from utils.config import (
    CLINVAR_VCF,
    GENOME_PATH,
    RESULTS_DIR,
)
from utils.gene_coordinates import GENES, get_gene
from utils.clinvar_parser import parse_clinvar_vcf, clinvar_to_variants
from utils.sequence_utils import (
    ScoringCheckpoint,
    Variant,
    generate_microhomology_deletions,
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
# ClinVar Indel Extraction
# =============================================================================

def extract_clinvar_indels(
    gene_symbol: str,
    min_stars: int = 1,
) -> Dict[str, List[Variant]]:
    """
    Extract ClinVar indels for a gene, separated by type.

    Returns: {
        "frameshift_plp": [...], "frameshift_blb": [...],
        "inframe_plp": [...], "inframe_blb": [...],
    }
    """
    result = {
        "frameshift_plp": [], "frameshift_blb": [],
        "inframe_plp": [], "inframe_blb": [],
    }

    for sig, sig_filter in [("plp", {"P/LP"}), ("blb", {"B/LB"})]:
        entries = parse_clinvar_vcf(
            str(CLINVAR_VCF), genes={gene_symbol},
            min_stars=min_stars, significance_filter=sig_filter,
        )
        variants = clinvar_to_variants(entries)

        for v in variants:
            if not v.is_indel:
                continue
            if v.is_frameshift:
                result["frameshift_{}".format(sig)].append(v)
            else:
                result["inframe_{}".format(sig)].append(v)

    return result


# =============================================================================
# Indel Analysis
# =============================================================================

def analyze_gene_indels(
    gene_symbol: str,
    window_size: int = 8192,
    min_stars: int = 1,
) -> Optional[dict]:
    """Analyze indel scoring results for a gene."""
    gene = get_gene(gene_symbol)
    log.info("\n--- {} ---".format(gene_symbol))

    # Load Evo2 results (includes indels if they were scored)
    ckpt_dir = RESULTS_DIR / gene_symbol / "w{}".format(window_size)
    checkpoint = ScoringCheckpoint(str(ckpt_dir))
    all_results = checkpoint.load_results()

    if not all_results:
        log.warning("  No Evo2 results")
        return None

    # Build delta lookup
    delta_lookup = {}
    for r in all_results:
        key = "{}:{}:{}>{}".format(r["chrom"], r["pos"], r["ref"], r["alt"])
        delta_lookup[key] = r["delta"]

    # Get ClinVar indels
    indels = extract_clinvar_indels(gene_symbol, min_stars=min_stars)

    result = {
        "gene": gene_symbol,
        "min_stars": min_stars,
        "counts": {},
        "metrics": {},
    }

    for category in ["frameshift", "inframe"]:
        plp_vars = indels["{}_plp".format(category)]
        blb_vars = indels["{}_blb".format(category)]

        # Match with Evo2 scores
        plp_deltas = []
        blb_deltas = []
        for v in plp_vars:
            if v.key() in delta_lookup:
                plp_deltas.append(delta_lookup[v.key()])
        for v in blb_vars:
            if v.key() in delta_lookup:
                blb_deltas.append(delta_lookup[v.key()])

        result["counts"][category] = {
            "clinvar_plp": len(plp_vars),
            "clinvar_blb": len(blb_vars),
            "scored_plp": len(plp_deltas),
            "scored_blb": len(blb_deltas),
        }

        log.info("  {}: ClinVar P/LP={}, B/LB={} | Scored P/LP={}, B/LB={}".format(
            category, len(plp_vars), len(blb_vars),
            len(plp_deltas), len(blb_deltas)))

        if plp_deltas:
            result["metrics"]["{}_plp".format(category)] = {
                "mean_delta": float(np.mean(plp_deltas)),
                "std_delta": float(np.std(plp_deltas)),
                "median_delta": float(np.median(plp_deltas)),
                "n": len(plp_deltas),
            }
            log.info("    {} P/LP: mean_delta={:.6f} (n={})".format(
                category, np.mean(plp_deltas), len(plp_deltas)))

        if blb_deltas:
            result["metrics"]["{}_blb".format(category)] = {
                "mean_delta": float(np.mean(blb_deltas)),
                "std_delta": float(np.std(blb_deltas)),
                "median_delta": float(np.median(blb_deltas)),
                "n": len(blb_deltas),
            }
            log.info("    {} B/LB: mean_delta={:.6f} (n={})".format(
                category, np.mean(blb_deltas), len(blb_deltas)))

        # Separation and AUROC if we have both classes
        if len(plp_deltas) >= 3 and len(blb_deltas) >= 3:
            separation = np.mean(blb_deltas) - np.mean(plp_deltas)
            result["metrics"]["{}_separation".format(category)] = float(separation)

            if scipy_stats:
                u_stat, u_p = scipy_stats.mannwhitneyu(
                    plp_deltas, blb_deltas, alternative="less"
                )
                result["metrics"]["{}_mannwhitney_p".format(category)] = float(u_p)
                log.info("    {} Mann-Whitney p={:.4e}".format(category, u_p))

            try:
                from utils.benchmarking import compute_auroc
                labels = np.array([1] * len(plp_deltas) + [0] * len(blb_deltas))
                scores = -np.array(plp_deltas + blb_deltas)
                auroc, _, _, _ = compute_auroc(labels, scores)
                result["metrics"]["{}_auroc".format(category)] = float(auroc)
                log.info("    {} AUROC={:.4f}".format(category, auroc))
            except ImportError:
                pass

    # Sanity check: frameshifts should have more negative deltas than in-frame
    fs_plp = result["metrics"].get("frameshift_plp", {})
    if_plp = result["metrics"].get("inframe_plp", {})
    if fs_plp and if_plp:
        log.info("  Frameshift vs in-frame P/LP mean delta: {:.6f} vs {:.6f}".format(
            fs_plp["mean_delta"], if_plp["mean_delta"]))
        result["sanity_check"] = {
            "frameshift_more_negative": fs_plp["mean_delta"] < if_plp["mean_delta"],
            "frameshift_mean_delta": fs_plp["mean_delta"],
            "inframe_mean_delta": if_plp["mean_delta"],
        }

    # Compare SNV vs indel distributions for scored variants
    snv_deltas = [r["delta"] for r in all_results
                  if len(r["ref"]) == 1 and len(r["alt"]) == 1]
    indel_deltas = [r["delta"] for r in all_results
                    if len(r["ref"]) != 1 or len(r["alt"]) != 1]

    if snv_deltas:
        result["snv_stats"] = {
            "n": len(snv_deltas),
            "mean_delta": float(np.mean(snv_deltas)),
            "std_delta": float(np.std(snv_deltas)),
        }
    if indel_deltas:
        result["indel_stats"] = {
            "n": len(indel_deltas),
            "mean_delta": float(np.mean(indel_deltas)),
            "std_delta": float(np.std(indel_deltas)),
        }

    return result


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Indel scoring analysis")
    parser.add_argument("--gene", type=str, default="all")
    parser.add_argument("--window-size", type=int, default=8192)
    parser.add_argument("--min-stars", type=int, default=1)
    args = parser.parse_args()

    all_symbols = [g.symbol for g in GENES.values()]
    genes = all_symbols if args.gene == "all" else [args.gene.upper()]

    log.info("Indel Analysis — genes: {}".format(genes))

    results = []
    for gene_sym in genes:
        result = analyze_gene_indels(gene_sym, args.window_size, args.min_stars)
        if result:
            results.append(result)

    if not results:
        log.error("No results computed")
        return

    # Summary
    log.info("\n" + "=" * 60)
    log.info("INDEL ANALYSIS SUMMARY")
    log.info("=" * 60)

    for r in results:
        gene = r["gene"]
        for cat in ["frameshift", "inframe"]:
            counts = r["counts"].get(cat, {})
            auroc = r["metrics"].get("{}_auroc".format(cat), None)
            auroc_str = "{:.4f}".format(auroc) if auroc is not None else "-"
            log.info("  {} {}: P/LP={} B/LB={} AUROC={}".format(
                gene, cat,
                counts.get("scored_plp", 0),
                counts.get("scored_blb", 0),
                auroc_str,
            ))

    # Save
    out_dir = RESULTS_DIR / "indel_analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(out_dir / "indel_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    # CSV for plotting
    import pandas as pd
    rows = []
    for r in results:
        for cat in ["frameshift", "inframe"]:
            for cls in ["plp", "blb"]:
                m = r["metrics"].get("{}_{}".format(cat, cls))
                if m:
                    rows.append({
                        "gene": r["gene"],
                        "category": cat,
                        "class": cls.upper(),
                        "mean_delta": m["mean_delta"],
                        "std_delta": m["std_delta"],
                        "n": m["n"],
                    })
    if rows:
        pd.DataFrame(rows).to_csv(out_dir / "indel_deltas.csv", index=False)

    log.info("\nResults saved to {}".format(out_dir))


if __name__ == "__main__":
    main()

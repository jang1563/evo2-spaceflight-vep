#!/usr/bin/env python3
"""
Radiation-characteristic mutation scoring (Phase 6).

Generates and scores radiation-characteristic mutations for each gene:
  - C>A and G>T SNVs (oxidative/proton/SBS18 signature via 8-oxoG)
  - Small deletions with microhomology (IR indel signature)

Compares functional impact of radiation-type mutations vs random mutations.
This is a two-step conditional analysis:
  "If space radiation produces these characteristic mutation types at this
   position, how functionally consequential would the resulting variant be?"

Produces: Figure 7 data

Usage:
  python 08_radiation_signatures.py --gene BRCA1 --window-size 8192
  python 08_radiation_signatures.py --gene all
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from utils.config import GENOME_PATH, RESULTS_DIR
from utils.gene_coordinates import GENES, get_gene
from utils.sequence_utils import (
    ScoringCheckpoint,
    Variant,
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
# Radiation Mutation Classification
# =============================================================================

# SBS18/oxidative damage signature: C>A and G>T transversions via 8-oxoguanine
# This is the dominant SNV signature from solar proton events and trapped protons.
# IMPORTANT: NOT C>T (that's UV/SBS7).

RADIATION_SNV_PAIRS = {
    # ref -> alt for radiation-characteristic transversions
    "C": "A",  # C>A via 8-oxoG on complementary strand
    "G": "T",  # G>T via 8-oxoG direct damage
}

# All possible SNV types for comparison
ALL_SNV_PAIRS = {
    "C": ["A", "G", "T"],
    "G": ["A", "C", "T"],
    "A": ["C", "G", "T"],
    "T": ["A", "C", "G"],
}


def classify_snv(ref: str, alt: str) -> str:
    """Classify an SNV by mutation type."""
    if ref in RADIATION_SNV_PAIRS and alt == RADIATION_SNV_PAIRS[ref]:
        return "radiation_oxidative"
    elif (ref == "C" and alt == "T") or (ref == "G" and alt == "A"):
        return "transition_CG_to_TA"  # UV-like, NOT radiation
    elif (ref == "A" and alt == "G") or (ref == "T" and alt == "C"):
        return "transition_AT_to_GC"
    else:
        return "other_transversion"


# =============================================================================
# Radiation Analysis
# =============================================================================

def analyze_radiation_mutations(
    gene_symbol: str,
    window_size: int = 8192,
) -> Optional[dict]:
    """
    Analyze radiation-characteristic mutations for a gene.

    Classifies all scored SNVs by mutation type and compares
    the functional impact (Evo2 delta) of radiation-type vs other mutations.
    """
    gene = get_gene(gene_symbol)
    log.info("\n--- {} ---".format(gene_symbol))

    # Load Evo2 results
    ckpt_dir = RESULTS_DIR / gene_symbol / "w{}".format(window_size)
    checkpoint = ScoringCheckpoint(str(ckpt_dir))
    results = checkpoint.load_results()

    if not results:
        log.warning("  No Evo2 results")
        return None

    # Classify each SNV
    by_class = {
        "radiation_oxidative": [],
        "transition_CG_to_TA": [],
        "transition_AT_to_GC": [],
        "other_transversion": [],
    }
    by_region = {
        "coding": {"radiation": [], "other": []},
        "intronic": {"radiation": [], "other": []},
        "promoter": {"radiation": [], "other": []},
        "utr5": {"radiation": [], "other": []},
        "utr3": {"radiation": [], "other": []},
    }

    for r in results:
        ref = r["ref"]
        alt = r["alt"]

        # SNVs only
        if len(ref) != 1 or len(alt) != 1:
            continue

        mut_class = classify_snv(ref, alt)
        delta = r["delta"]

        if mut_class in by_class:
            by_class[mut_class].append(delta)

        # Region-specific analysis
        region = r.get("region_type", "")
        is_radiation = mut_class == "radiation_oxidative"

        for reg_name in by_region:
            if region == reg_name or (region == "" and reg_name == "coding"):
                key = "radiation" if is_radiation else "other"
                by_region[reg_name][key].append(delta)

    # Indel analysis (microhomology deletions)
    indel_deltas = {
        "microhomology_del": [],
        "other_indel": [],
    }
    for r in results:
        ref = r["ref"]
        alt = r["alt"]
        if len(ref) == 1 and len(alt) == 1:
            continue  # Skip SNVs
        region_type = r.get("region_type", "")
        if "radiation_del_mh" in region_type:
            indel_deltas["microhomology_del"].append(r["delta"])
        elif len(ref) != len(alt):
            indel_deltas["other_indel"].append(r["delta"])

    # Compute statistics
    result = {
        "gene": gene_symbol,
        "snv_classes": {},
        "region_analysis": {},
        "indel_analysis": {},
    }

    for cls_name, deltas in by_class.items():
        if not deltas:
            continue
        arr = np.array(deltas)
        result["snv_classes"][cls_name] = {
            "n": len(deltas),
            "mean_delta": float(np.mean(arr)),
            "std_delta": float(np.std(arr)),
            "median_delta": float(np.median(arr)),
            "q10": float(np.percentile(arr, 10)),
            "q90": float(np.percentile(arr, 90)),
        }
        log.info("  {}: n={}, mean_delta={:.6f}".format(
            cls_name, len(deltas), np.mean(arr)))

    # Compare radiation vs non-radiation
    rad_deltas = by_class.get("radiation_oxidative", [])
    nonrad_deltas = (
        by_class.get("transition_CG_to_TA", [])
        + by_class.get("transition_AT_to_GC", [])
        + by_class.get("other_transversion", [])
    )

    if rad_deltas and nonrad_deltas and scipy_stats:
        rad_arr = np.array(rad_deltas)
        nonrad_arr = np.array(nonrad_deltas)
        u_stat, u_p = scipy_stats.mannwhitneyu(rad_arr, nonrad_arr, alternative="two-sided")
        ks_stat, ks_p = scipy_stats.ks_2samp(rad_arr, nonrad_arr)

        result["radiation_vs_other"] = {
            "radiation_n": len(rad_deltas),
            "other_n": len(nonrad_deltas),
            "radiation_mean_delta": float(np.mean(rad_arr)),
            "other_mean_delta": float(np.mean(nonrad_arr)),
            "mannwhitney_p": float(u_p),
            "ks_stat": float(ks_stat),
            "ks_p": float(ks_p),
        }
        log.info("  Radiation vs other: mean_delta {:.6f} vs {:.6f}, "
                 "MWU p={:.4e}".format(
                     np.mean(rad_arr), np.mean(nonrad_arr), u_p))

    # Region-specific comparison
    for region, data in by_region.items():
        if data["radiation"] and data["other"]:
            rad_arr = np.array(data["radiation"])
            other_arr = np.array(data["other"])
            entry = {
                "radiation_n": len(data["radiation"]),
                "other_n": len(data["other"]),
                "radiation_mean": float(np.mean(rad_arr)),
                "other_mean": float(np.mean(other_arr)),
            }
            if scipy_stats and len(data["radiation"]) >= 5:
                _, p = scipy_stats.mannwhitneyu(rad_arr, other_arr, alternative="two-sided")
                entry["mannwhitney_p"] = float(p)
            result["region_analysis"][region] = entry

    # Indel results
    for indel_type, deltas in indel_deltas.items():
        if deltas:
            arr = np.array(deltas)
            result["indel_analysis"][indel_type] = {
                "n": len(deltas),
                "mean_delta": float(np.mean(arr)),
                "std_delta": float(np.std(arr)),
                "median_delta": float(np.median(arr)),
            }

    # Identify radiation-sensitive functional hotspots
    # Positions where radiation-type mutations have the most negative deltas
    rad_variants = []
    for r in results:
        if len(r["ref"]) == 1 and len(r["alt"]) == 1:
            cls = classify_snv(r["ref"], r["alt"])
            if cls == "radiation_oxidative":
                rad_variants.append(r)

    if rad_variants:
        # Sort by delta (most negative = most damaging)
        rad_variants.sort(key=lambda x: x["delta"])
        top_n = min(20, len(rad_variants))
        result["top_radiation_hotspots"] = [
            {
                "variant_key": "{}:{}:{}>{}".format(
                    r["chrom"], r["pos"], r["ref"], r["alt"]),
                "delta": float(r["delta"]),
                "region_type": r.get("region_type", ""),
                "pos": r["pos"],
            }
            for r in rad_variants[:top_n]
        ]

    return result


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Radiation signature mutation scoring"
    )
    parser.add_argument("--gene", type=str, default="all")
    parser.add_argument("--window-size", type=int, default=8192)
    args = parser.parse_args()

    all_symbols = [g.symbol for g in GENES.values()]
    genes = all_symbols if args.gene == "all" else [args.gene.upper()]

    log.info("Radiation Signature Analysis — genes: {}".format(genes))

    results = []
    for gene_sym in genes:
        result = analyze_radiation_mutations(gene_sym, args.window_size)
        if result:
            results.append(result)

    if not results:
        log.error("No results computed")
        return

    # Aggregate summary
    log.info("\n" + "=" * 60)
    log.info("RADIATION SIGNATURE SUMMARY")
    log.info("=" * 60)

    log.info("{:<8} {:>8} {:>12} {:>12} {:>10}".format(
        "Gene", "N_rad", "Rad_mean_d", "Other_mean_d", "MWU_p"))
    log.info("-" * 52)

    for r in results:
        comp = r.get("radiation_vs_other", {})
        if comp:
            log.info("{:<8} {:>8} {:>12.6f} {:>12.6f} {:>10.2e}".format(
                r["gene"],
                comp.get("radiation_n", 0),
                comp.get("radiation_mean_delta", 0),
                comp.get("other_mean_delta", 0),
                comp.get("mannwhitney_p", 1),
            ))

    # Save results
    out_dir = RESULTS_DIR / "radiation_signatures"
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(out_dir / "radiation_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    # Per-variant data for violin plots (Figure 7C)
    import pandas as pd
    rows = []
    for r in results:
        for cls_name, stats in r["snv_classes"].items():
            rows.append({
                "gene": r["gene"],
                "mutation_class": cls_name,
                "mean_delta": stats["mean_delta"],
                "std_delta": stats["std_delta"],
                "n": stats["n"],
            })
    if rows:
        pd.DataFrame(rows).to_csv(
            out_dir / "radiation_class_summary.csv", index=False
        )

    # Top hotspots across all genes
    all_hotspots = []
    for r in results:
        for hs in r.get("top_radiation_hotspots", []):
            hs["gene"] = r["gene"]
            all_hotspots.append(hs)

    if all_hotspots:
        pd.DataFrame(all_hotspots).to_csv(
            out_dir / "radiation_hotspots.csv", index=False
        )

    log.info("\nResults saved to {}".format(out_dir))


if __name__ == "__main__":
    main()

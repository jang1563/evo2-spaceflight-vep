#!/usr/bin/env python3
"""
Score published astronaut mutations with Evo2 (Phase 7).

Scores specific variants identified in astronauts (Brojakowska 2022,
I4 multi-omics 2024) and contextualizes them against ClinVar and DMS data.

Produces: Figure 8 data

Usage:
  python 11_astronaut_variants.py --window-size 8192
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from utils.config import RESULTS_DIR
from utils.gene_coordinates import GENES, get_gene
from utils.sequence_utils import ScoringCheckpoint, Variant

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)


# =============================================================================
# Known Astronaut Variants (Brojakowska 2022, Commun Biol)
# =============================================================================

# These will be populated from supplementary data.
# For now, define the structure and known hotspot variants.
# Format: (gene, chrom, pos_0based, ref, alt, note)

KNOWN_ASTRONAUT_MUTATIONS = [
    # TP53 — 7 variants in 14 astronauts (most mutated CH gene)
    # Positions verified via CDS map from NM_000546.6 (GRCh38)
    ("TP53", "chr17", 7675087, "C", "T", "R175H — most common TP53 mutation"),
    ("TP53", "chr17", 7674220, "G", "A", "R248W — TP53 DNA-binding hotspot"),
    ("TP53", "chr17", 7673801, "C", "T", "R273H — TP53 DNA-contact residue"),
    ("TP53", "chr17", 7673775, "G", "A", "R282W — TP53 structural mutation"),

    # DNMT3A — 6 variants in 14 astronauts (2nd most mutated)
    # Positions verified via CDS map from NM_175629.2 (GRCh38)
    ("DNMT3A", "chr2", 25234372, "C", "T", "R882H — canonical DNMT3A CH hotspot"),
    ("DNMT3A", "chr2", 25234373, "G", "A", "R882C — 2nd most common DNMT3A CH"),

    # TERT promoter — not from astronaut CH, but relevant for
    # radiation-induced somatic mutations in long-duration spaceflight
    # TERT is minus strand. C228T/C250T are named by minus-strand position.
    # On the forward strand (which Evo2 scores), both are G>A substitutions.
    # GRCh38 coordinates verified against reference genome and MPRA mapping.
    ("TERT", "chr5", 1295112, "G", "A", "C228T — TERT promoter hotspot (creates ETS site)"),
    ("TERT", "chr5", 1295134, "G", "A", "C250T — TERT promoter hotspot (creates ETS site)"),
]


def score_astronaut_variants(
    window_size: int = 8192,
) -> List[dict]:
    """
    Score known astronaut variants and contextualize them.

    For each variant:
    1. Look up Evo2 score from existing results
    2. Compare against gene-wide score distribution
    3. Rank among ClinVar variants
    """
    results = []

    for gene_sym, chrom, pos, ref, alt, note in KNOWN_ASTRONAUT_MUTATIONS:
        variant_key = "{}:{}:{}>{}".format(chrom, pos, ref, alt)
        log.info("  {} — {}".format(variant_key, note))

        # Load Evo2 results for this gene
        ckpt_dir = RESULTS_DIR / gene_sym / "w{}".format(window_size)
        checkpoint = ScoringCheckpoint(str(ckpt_dir))
        all_results = checkpoint.load_results()

        if not all_results:
            log.warning("    No Evo2 results for {}".format(gene_sym))
            results.append({
                "gene": gene_sym,
                "variant_key": variant_key,
                "note": note,
                "evo2_delta": None,
                "status": "not_scored",
            })
            continue

        # Build delta lookup
        delta_lookup = {}
        all_deltas = []
        for r in all_results:
            key = "{}:{}:{}>{}".format(r["chrom"], r["pos"], r["ref"], r["alt"])
            delta_lookup[key] = r["delta"]
            all_deltas.append(r["delta"])

        all_deltas = np.array(all_deltas)

        entry = {
            "gene": gene_sym,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "variant_key": variant_key,
            "note": note,
        }

        if variant_key in delta_lookup:
            delta = delta_lookup[variant_key]
            entry["evo2_delta"] = float(delta)
            entry["evo2_pathogenicity"] = float(-delta)
            entry["status"] = "scored"

            # Percentile rank among all gene variants
            percentile = float(np.mean(all_deltas >= delta) * 100)
            entry["percentile_rank"] = percentile
            entry["total_variants_scored"] = len(all_deltas)

            # Compare to ClinVar classes
            plp_deltas = [r["delta"] for r in all_results
                          if r.get("clinvar_class") == "P/LP"]
            blb_deltas = [r["delta"] for r in all_results
                          if r.get("clinvar_class") == "B/LB"]

            if plp_deltas:
                entry["plp_mean_delta"] = float(np.mean(plp_deltas))
                entry["more_damaging_than_avg_plp"] = delta < np.mean(plp_deltas)
            if blb_deltas:
                entry["blb_mean_delta"] = float(np.mean(blb_deltas))

            log.info("    delta={:.6f}, percentile={:.1f}%".format(delta, percentile))
        else:
            entry["evo2_delta"] = None
            entry["status"] = "not_in_results"
            log.info("    Not found in scored variants")

        results.append(entry)

    return results


def contextualize_with_clinvar(
    astronaut_results: List[dict],
    window_size: int = 8192,
) -> List[dict]:
    """
    Add ClinVar context to astronaut variant results.

    For each astronaut variant, find its ClinVar classification (if any).
    """
    from utils.clinvar_parser import parse_clinvar_vcf
    from utils.config import CLINVAR_VCF

    # Build ClinVar lookup
    clinvar_lookup = {}
    genes_needed = set(r["gene"] for r in astronaut_results)

    for gene_sym in genes_needed:
        entries = parse_clinvar_vcf(
            str(CLINVAR_VCF), genes={gene_sym}, min_stars=0,
        )
        for e in entries:
            key = "{}:{}:{}>{}".format(e.chrom, e.pos, e.ref, e.alt)
            clinvar_lookup[key] = {
                "significance": e.significance,
                "stars": e.stars,
                "clinvar_id": e.clinvar_id,
            }

    for result in astronaut_results:
        key = result["variant_key"]
        if key in clinvar_lookup:
            cv = clinvar_lookup[key]
            result["clinvar_significance"] = cv["significance"]
            result["clinvar_stars"] = cv["stars"]
            result["clinvar_id"] = cv["clinvar_id"]
        else:
            result["clinvar_significance"] = "Not in ClinVar"

    return astronaut_results


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Score astronaut variants")
    parser.add_argument("--window-size", type=int, default=8192)
    args = parser.parse_args()

    log.info("=" * 60)
    log.info("Astronaut Variant Scoring")
    log.info("=" * 60)

    # Score variants
    results = score_astronaut_variants(args.window_size)

    # Add ClinVar context
    try:
        results = contextualize_with_clinvar(results, args.window_size)
    except Exception as e:
        log.warning("ClinVar contextualization failed: {}".format(e))

    # Summary
    log.info("\n" + "=" * 60)
    log.info("ASTRONAUT VARIANT SUMMARY")
    log.info("=" * 60)

    log.info("{:<8} {:<30} {:>10} {:>8} {:<12}".format(
        "Gene", "Variant", "Delta", "Pctile", "ClinVar"))
    log.info("-" * 70)

    for r in results:
        delta_str = "{:.6f}".format(r["evo2_delta"]) if r.get("evo2_delta") is not None else "-"
        pctile_str = "{:.1f}%".format(r.get("percentile_rank", 0)) if r.get("evo2_delta") is not None else "-"
        clinvar = r.get("clinvar_significance", "?")
        log.info("{:<8} {:<30} {:>10} {:>8} {:<12}".format(
            r["gene"], r["note"][:30], delta_str, pctile_str, clinvar))

    # Save
    out_dir = RESULTS_DIR / "astronaut_variants"
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(out_dir / "astronaut_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    import pandas as pd
    pd.DataFrame(results).to_csv(out_dir / "astronaut_variants.csv", index=False)

    log.info("\nResults saved to {}".format(out_dir))

    # Scored count
    scored = sum(1 for r in results if r.get("evo2_delta") is not None)
    log.info("{}/{} variants scored".format(scored, len(results)))


if __name__ == "__main__":
    main()

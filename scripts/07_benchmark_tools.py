#!/usr/bin/env python3
"""
Multi-tool benchmark comparison (Phase 3).

Compare Evo2 against AlphaMissense, CADD, REVEL, and SpliceAI
on ClinVar variants across all genes.

Computes:
- Per-gene AUROC/AUPRC for each tool
- Pairwise Spearman correlations (orthogonality matrix)
- Tool coverage statistics
- Concordance/discordance analysis

Produces: Figure 4 (tool comparison bars), Supp Fig S1 (orthogonality matrix)

Usage:
  python 07_benchmark_tools.py --gene all
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
    compute_auroc,
    compute_auprc,
    compute_orthogonality_matrix,
    interpret_orthogonality,
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


def load_all_tool_scores(gene_symbol: str, variant_keys: List[str]) -> Dict[str, Dict[str, float]]:
    """
    Load scores from all benchmark tools for a set of variants.

    Returns: {tool_name: {variant_key: score}}
    """
    gene = get_gene(gene_symbol)
    tool_scores = {}

    # AlphaMissense (via tabix)
    try:
        import pysam
        am_path = str(ALPHAMISSENSE_TSV)
        if Path(am_path).exists():
            scores = {}
            tbx = pysam.TabixFile(am_path)
            chrom_q = gene.chrom if gene.chrom.startswith("chr") else "chr" + gene.chrom
            for row in tbx.fetch(chrom_q, gene.start, gene.end):
                fields = row.split("\t")
                chrom = fields[0] if fields[0].startswith("chr") else "chr" + fields[0]
                pos = int(fields[1]) - 1
                key = "{}:{}:{}>{}".format(chrom, pos, fields[2], fields[3])
                if key in variant_keys:
                    scores[key] = float(fields[8])
            tbx.close()
            if scores:
                tool_scores["AlphaMissense"] = scores
                log.info("    AlphaMissense: {} variants".format(len(scores)))
    except Exception as e:
        log.warning("    AlphaMissense load failed: {}".format(e))

    # CADD v1.7 (via tabix)
    try:
        import pysam
        cadd_path = str(SHARED_DATA_DIR / "benchmarks" / "cadd" / "whole_genome_SNVs.tsv.gz")
        if Path(cadd_path).exists():
            scores = {}
            tbx = pysam.TabixFile(cadd_path)
            chrom_q = gene.chrom.replace("chr", "")
            for row in tbx.fetch(chrom_q, gene.start, gene.end):
                fields = row.split("\t")
                chrom = "chr" + fields[0]
                pos = int(fields[1]) - 1
                key = "{}:{}:{}>{}".format(chrom, pos, fields[2], fields[3])
                if key in variant_keys:
                    scores[key] = float(fields[5])  # PHRED score
            tbx.close()
            if scores:
                tool_scores["CADD"] = scores
                log.info("    CADD: {} variants".format(len(scores)))
    except Exception as e:
        log.warning("    CADD load failed: {}".format(e))

    # REVEL (via tabix on processed BED)
    try:
        import pysam
        revel_path = str(SHARED_DATA_DIR / "benchmarks" / "revel" / "revel_genes.bed.gz")
        if Path(revel_path).exists():
            scores = {}
            tbx = pysam.TabixFile(revel_path)
            chrom_q = gene.chrom if gene.chrom.startswith("chr") else "chr" + gene.chrom
            for row in tbx.fetch(chrom_q, gene.start, gene.end):
                fields = row.split("\t")
                # BED format: chrom, start, end, ref, alt, aaref, aaalt, REVEL
                chrom = fields[0]
                pos = int(fields[1])  # 0-based BED start
                ref = fields[3]
                alt = fields[4]
                revel_score = float(fields[7])
                key = "{}:{}:{}>{}".format(chrom, pos, ref, alt)
                if key in variant_keys:
                    scores[key] = revel_score
            tbx.close()
            if scores:
                tool_scores["REVEL"] = scores
                log.info("    REVEL: {} variants".format(len(scores)))
    except Exception as e:
        log.warning("    REVEL load failed: {}".format(e))

    # SpliceAI (from extracted per-gene TSV)
    try:
        import gzip
        spliceai_path = str(
            SHARED_DATA_DIR / "benchmarks" / "spliceai"
            / "{}_spliceai_scores.tsv.gz".format(gene_symbol)
        )
        if Path(spliceai_path).exists():
            scores = {}
            with gzip.open(spliceai_path, "rt") as fh:
                header = fh.readline()  # skip header
                for line in fh:
                    fields = line.rstrip("\n").split("\t")
                    chrom = fields[0]
                    pos = int(fields[1]) - 1  # VCF 1-based → 0-based
                    ref = fields[2]
                    alt = fields[3]
                    # Max delta score across the 4 splice change types
                    ds_max = max(
                        float(fields[5]), float(fields[6]),
                        float(fields[7]), float(fields[8]),
                    )
                    key = "{}:{}:{}>{}".format(chrom, pos, ref, alt)
                    if key in variant_keys:
                        scores[key] = ds_max
            if scores:
                tool_scores["SpliceAI"] = scores
                log.info("    SpliceAI: {} variants".format(len(scores)))
    except Exception as e:
        log.warning("    SpliceAI load failed: {}".format(e))

    # ncER (10bp bins, hg38, from processed BED)
    try:
        import pysam
        ncer_path = str(
            SHARED_DATA_DIR / "noncoding" / "ncer"
            / "ncER_hg38_gene_regions.bed.gz"
        )
        if Path(ncer_path).exists():
            scores = {}
            tbx = pysam.TabixFile(ncer_path)
            chrom_q = gene.chrom if gene.chrom.startswith("chr") else "chr" + gene.chrom
            # Build position → ncER lookup from 10bp bins
            bin_scores = {}
            for row in tbx.fetch(chrom_q, gene.start, gene.end):
                fields = row.split("\t")
                bin_start = int(fields[1])
                bin_end = int(fields[2])
                ncer_score = float(fields[3])
                for p in range(bin_start, bin_end):
                    bin_scores[p] = ncer_score
            tbx.close()
            # Assign to variants
            for key in variant_keys:
                parts = key.split(":")
                pos = int(parts[1])
                if pos in bin_scores:
                    scores[key] = bin_scores[pos]
            if scores:
                tool_scores["ncER"] = scores
                log.info("    ncER: {} variants".format(len(scores)))
    except Exception as e:
        log.warning("    ncER load failed: {}".format(e))

    return tool_scores


def compute_gene_comparison(
    gene_symbol: str,
    window_size: int = 8192,
    min_stars: int = 2,
) -> Optional[dict]:
    """Compute full tool comparison for one gene."""
    gene = get_gene(gene_symbol)
    log.info("\n--- {} ---".format(gene_symbol))

    # Load Evo2 results
    ckpt_dir = RESULTS_DIR / gene_symbol / "w{}".format(window_size)
    checkpoint = ScoringCheckpoint(str(ckpt_dir))
    evo2_results = checkpoint.load_results()

    if not evo2_results:
        log.warning("  No Evo2 results")
        return None

    # Build Evo2 delta lookup
    evo2_deltas = {}
    for r in evo2_results:
        key = "{}:{}:{}>{}".format(r["chrom"], r["pos"], r["ref"], r["alt"])
        evo2_deltas[key] = r["delta"]

    # Get ClinVar variants ≥min_stars
    entries_plp = parse_clinvar_vcf(
        str(CLINVAR_VCF), genes={gene_symbol},
        min_stars=min_stars, significance_filter={"P/LP"},
    )
    entries_blb = parse_clinvar_vcf(
        str(CLINVAR_VCF), genes={gene_symbol},
        min_stars=min_stars, significance_filter={"B/LB"},
    )

    # SNVs only, matched to Evo2
    plp_keys = []
    blb_keys = []
    for e in clinvar_to_variants(entries_plp):
        if e.is_snv and e.key() in evo2_deltas:
            plp_keys.append(e.key())
    for e in clinvar_to_variants(entries_blb):
        if e.is_snv and e.key() in evo2_deltas:
            blb_keys.append(e.key())

    all_keys = plp_keys + blb_keys
    all_labels = [1] * len(plp_keys) + [0] * len(blb_keys)

    if len(plp_keys) < 5 or len(blb_keys) < 5:
        log.warning("  Too few ClinVar variants (P/LP={}, B/LB={})".format(
            len(plp_keys), len(blb_keys)))
        return None

    log.info("  ClinVar ≥{}★: {} P/LP, {} B/LB".format(
        min_stars, len(plp_keys), len(blb_keys)))

    # Load all tool scores
    variant_key_set = set(all_keys)
    tool_scores = load_all_tool_scores(gene_symbol, variant_key_set)

    # Add Evo2
    evo2_scores = {}
    for k in all_keys:
        evo2_scores[k] = -evo2_deltas[k]  # Negate: higher = more pathogenic
    tool_scores["Evo2"] = evo2_scores

    # Compute AUROC for each tool
    tool_metrics = {}
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

        auroc, _, _, _ = compute_auroc(t_labels, t_scores)
        auprc, _, _, _ = compute_auprc(t_labels, t_scores)

        tool_metrics[tool_name] = {
            "auroc": float(auroc),
            "auprc": float(auprc),
            "n_scored": len(t_labels),
            "n_plp": int(sum(t_labels)),
            "n_blb": int(sum(1 - t_labels)),
        }
        log.info("  {}: AUROC={:.4f}, AUPRC={:.4f} (n={})".format(
            tool_name, auroc, auprc, len(t_labels)))

    # Pairwise orthogonality on shared variants
    ortho = {}
    if scipy_stats and len(tool_scores) >= 2:
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
                pair = "{}_vs_{}".format(t1, t2)
                ortho[pair] = {
                    "rho": float(rho),
                    "p": float(p),
                    "n": len(shared),
                    "interpretation": interpret_orthogonality(rho),
                }

    return {
        "gene": gene_symbol,
        "min_stars": min_stars,
        "n_plp": len(plp_keys),
        "n_blb": len(blb_keys),
        "tool_metrics": tool_metrics,
        "orthogonality": ortho,
    }


def main():
    parser = argparse.ArgumentParser(description="Multi-tool benchmark comparison")
    parser.add_argument("--gene", type=str, default="all")
    parser.add_argument("--window-size", type=int, default=8192)
    parser.add_argument("--min-stars", type=int, default=2)
    args = parser.parse_args()

    all_symbols = [g.symbol for g in GENES.values()]
    genes = all_symbols if args.gene == "all" else [args.gene.upper()]

    log.info("Tool Benchmark — genes: {}".format(genes))

    results = []
    for gene_sym in genes:
        result = compute_gene_comparison(gene_sym, args.window_size, args.min_stars)
        if result:
            results.append(result)

    if not results:
        log.error("No results computed")
        return

    # Aggregate orthogonality across genes
    log.info("\n{}".format("=" * 80))
    log.info("AGGREGATE RESULTS")
    log.info("=" * 80)

    # Per-gene AUROC table
    all_tools = set()
    for r in results:
        all_tools.update(r["tool_metrics"].keys())
    all_tools = sorted(all_tools)

    header = "{:<8}".format("Gene")
    for t in all_tools:
        header += " {:>12}".format(t[:12])
    log.info(header)
    log.info("-" * (8 + 13 * len(all_tools)))

    for r in results:
        line = "{:<8}".format(r["gene"])
        for t in all_tools:
            m = r["tool_metrics"].get(t)
            if m:
                line += " {:>12.4f}".format(m["auroc"])
            else:
                line += " {:>12}".format("-")
        log.info(line)

    # Mean AUROC per tool
    log.info("")
    log.info("Mean AUROC across genes:")
    for t in all_tools:
        aurocs = [r["tool_metrics"][t]["auroc"]
                  for r in results if t in r["tool_metrics"]]
        if aurocs:
            log.info("  {}: {:.4f} ± {:.4f} (n={} genes)".format(
                t, np.mean(aurocs), np.std(aurocs), len(aurocs)))

    # Save
    out_dir = RESULTS_DIR / "benchmark_comparison"
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(out_dir / "benchmark_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    # Save as CSV for plotting
    rows = []
    for r in results:
        for tool, m in r["tool_metrics"].items():
            rows.append({
                "gene": r["gene"],
                "tool": tool,
                "auroc": m["auroc"],
                "auprc": m["auprc"],
                "n_scored": m["n_scored"],
            })
    pd.DataFrame(rows).to_csv(out_dir / "tool_aurocs.csv", index=False)

    log.info("\nSaved to {}".format(out_dir))


if __name__ == "__main__":
    main()

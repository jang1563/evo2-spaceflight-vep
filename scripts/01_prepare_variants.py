#!/usr/bin/env python3
"""
Prepare variant sets for Evo2 scoring.

For each gene, generates:
1. ClinVar variants (P/LP, B/LB, VUS) with star stratification
2. Saturation mutagenesis (all possible coding SNVs)
3. Non-coding SNVs (promoter + intronic, subsampled if too large)
4. DMS-matched variants (for control genes)
5. Radiation-characteristic mutations (C>A, microhomology deletions)

Output: per-gene JSONL files in results/<GENE>/variants.jsonl

Usage:
  python 01_prepare_variants.py --gene BRCA1
  python 01_prepare_variants.py --gene all
"""

import argparse
import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from utils.config import CLINVAR_VCF, GENOME_PATH, RESULTS_DIR
from utils.gene_coordinates import GENES, get_gene, CONTROL_GENES
from utils.sequence_utils import (
    GenomeAccessor,
    Variant,
    annotate_variant_region,
    generate_all_snvs,
    generate_radiation_snvs,
    generate_microhomology_deletions,
)
from utils.clinvar_parser import parse_clinvar_vcf, clinvar_to_variants

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)


def save_variants(variants: list, filepath: str):
    """Save variants to JSONL file."""
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, 'w') as f:
        for v in variants:
            data = {
                "chrom": v.chrom,
                "pos": v.pos,
                "ref": v.ref,
                "alt": v.alt,
                "variant_id": v.variant_id,
                "gene": v.gene,
                "region_type": v.region_type,
                "clinvar_class": v.clinvar_class,
                "clinvar_stars": v.clinvar_stars,
                "dms_score": v.dms_score,
            }
            f.write(json.dumps(data) + "\n")


def prepare_gene(gene_symbol: str, genome: GenomeAccessor):
    """Prepare all variant sets for a single gene."""
    gene = get_gene(gene_symbol)
    out_dir = RESULTS_DIR / gene_symbol
    out_dir.mkdir(parents=True, exist_ok=True)

    all_variants = []
    variant_keys = set()

    def add_unique(variants):
        """Add variants, deduplicating by key."""
        added = 0
        for v in variants:
            k = v.key()
            if k not in variant_keys:
                variant_keys.add(k)
                all_variants.append(v)
                added += 1
        return added

    # 1. ClinVar variants
    log.info(f"  Fetching ClinVar variants for {gene_symbol}...")
    clinvar_entries = parse_clinvar_vcf(
        str(CLINVAR_VCF),
        genes={gene_symbol},
    )
    clinvar_vars = clinvar_to_variants(clinvar_entries)
    for v in clinvar_vars:
        v.gene = gene_symbol
        v.region_type = annotate_variant_region(v, gene)
    n_cv = add_unique(clinvar_vars)
    log.info(f"    ClinVar: {n_cv} variants "
             f"(P/LP: {sum(1 for v in clinvar_vars if v.clinvar_class == 'P/LP')}, "
             f"B/LB: {sum(1 for v in clinvar_vars if v.clinvar_class == 'B/LB')}, "
             f"VUS: {sum(1 for v in clinvar_vars if v.clinvar_class == 'VUS')})")

    # 2. Saturation coding SNVs
    log.info(f"  Generating saturation coding SNVs for {gene_symbol}...")
    coding_snvs = []
    for exon in gene.coding_exons:
        ex_start = max(exon.start, gene.cds_start)
        ex_end = min(exon.end, gene.cds_end)
        exon_snvs = generate_all_snvs(genome, gene.chrom, ex_start, ex_end, gene_symbol)
        for v in exon_snvs:
            v.region_type = "coding"
        coding_snvs.extend(exon_snvs)
    n_coding = add_unique(coding_snvs)
    log.info(f"    Coding SNVs: {n_coding}")

    # 3. Promoter SNVs
    log.info(f"  Generating promoter SNVs for {gene_symbol}...")
    prom_chrom, prom_start, prom_end = gene.get_promoter_region(upstream=2000, downstream=500)
    prom_snvs = generate_all_snvs(genome, prom_chrom, prom_start, prom_end, gene_symbol)
    for v in prom_snvs:
        v.region_type = "promoter"
    n_prom = add_unique(prom_snvs)
    log.info(f"    Promoter SNVs: {n_prom}")

    # 4. Radiation-characteristic mutations (coding + promoter)
    log.info(f"  Generating radiation-characteristic mutations for {gene_symbol}...")
    # C>A / G>T in coding regions
    rad_snvs = []
    for exon in gene.coding_exons:
        ex_start = max(exon.start, gene.cds_start)
        ex_end = min(exon.end, gene.cds_end)
        rad_snvs.extend(generate_radiation_snvs(genome, gene.chrom, ex_start, ex_end, gene_symbol))
    # C>A / G>T in promoter
    rad_snvs.extend(generate_radiation_snvs(genome, prom_chrom, prom_start, prom_end, gene_symbol))
    n_rad_snv = add_unique(rad_snvs)

    # Microhomology deletions in coding
    mh_dels = []
    for exon in gene.coding_exons:
        ex_start = max(exon.start, gene.cds_start)
        ex_end = min(exon.end, gene.cds_end)
        mh_dels.extend(generate_microhomology_deletions(
            genome, gene.chrom, ex_start, ex_end, gene_symbol
        ))
    n_mh = add_unique(mh_dels)
    log.info(f"    Radiation SNVs: {n_rad_snv}, MH deletions: {n_mh}")

    # 5. Summary
    log.info(f"  Total unique variants for {gene_symbol}: {len(all_variants)}")

    # Save
    out_file = out_dir / "variants.jsonl"
    save_variants(all_variants, str(out_file))
    log.info(f"  Saved to {out_file}")

    # Also save summary
    summary = {
        "gene": gene_symbol,
        "total_variants": len(all_variants),
        "clinvar": n_cv,
        "coding_snvs": n_coding,
        "promoter_snvs": n_prom,
        "radiation_snvs": n_rad_snv,
        "microhomology_dels": n_mh,
    }
    with open(out_dir / "variant_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    return summary


def main():
    parser = argparse.ArgumentParser(description="Prepare variants for Evo2 scoring")
    parser.add_argument("--gene", type=str, default="all",
                        help="Gene symbol or 'all'")
    args = parser.parse_args()

    genes = list(GENES.keys()) if args.gene == "all" else [args.gene.upper()]

    log.info(f"Preparing variants for {len(genes)} gene(s): {genes}")

    genome = GenomeAccessor(str(GENOME_PATH))

    summaries = []
    for gene_sym in genes:
        log.info(f"\n{'='*60}")
        log.info(f"Processing {gene_sym}")
        log.info(f"{'='*60}")
        summary = prepare_gene(gene_sym, genome)
        summaries.append(summary)

    genome.close()

    # Print overall summary
    log.info(f"\n{'='*60}")
    log.info("Overall Summary")
    log.info(f"{'='*60}")
    total = 0
    for s in summaries:
        log.info(f"  {s['gene']}: {s['total_variants']:,} variants")
        total += s["total_variants"]
    log.info(f"  TOTAL: {total:,} variants")


if __name__ == "__main__":
    main()

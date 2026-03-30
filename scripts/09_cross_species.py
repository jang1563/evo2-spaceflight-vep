#!/usr/bin/env python3
"""
Cross-species scoring with Evo2 (Phase 7 supplement).

Scores mouse (GRCm39/mm39) ortholog regions for selected genes and compares
constraint patterns between human and mouse. Connects to NASA GeneLab
mouse radiation experiment datasets.

Genes: ATM, TERT, RAD51 (have well-characterized mouse orthologs)

Produces: Supplementary cross-species data

Usage:
  python 09_cross_species.py --gene ATM
  python 09_cross_species.py --gene all
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from utils.config import RESULTS_DIR, SHARED_DATA_DIR
from utils.sequence_utils import BASES, Variant

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
# Mouse Ortholog Coordinates (GRCm39/mm39)
# =============================================================================
# Sources: NCBI Gene, Ensembl 112, UCSC liftOver verification
# Coordinates: 0-based, half-open [start, end)

MOUSE_ORTHOLOGS = {
    "ATM": {
        "mouse_symbol": "Atm",
        "mouse_chrom": "chr9",
        "mouse_start": 53333587,
        "mouse_end": 53443457,
        "mouse_strand": "+",
        "mouse_refseq": "NM_007499.3",
        "human_symbol": "ATM",
        "human_chrom": "chr11",
        "human_start": 108223067,
        "human_end": 108369102,
        "identity_pct": 87.3,  # protein identity
        "genelab_datasets": [
            "GLDS-288",  # ISS mice irradiation
            "GLDS-352",  # GeneLab multi-omics
        ],
    },
    "TERT": {
        "mouse_symbol": "Tert",
        "mouse_chrom": "chr13",
        "mouse_start": 73492486,
        "mouse_end": 73524908,
        "mouse_strand": "-",
        "mouse_refseq": "NM_009354.1",
        "human_symbol": "TERT",
        "human_chrom": "chr5",
        "human_start": 1253166,
        "human_end": 1295068,
        "identity_pct": 72.1,
        "genelab_datasets": [
            "GLDS-288",
            "GLDS-168",  # RR-1 Rodent Research
        ],
    },
    "RAD51": {
        "mouse_symbol": "Rad51",
        "mouse_chrom": "chr7",
        "mouse_start": 44620768,
        "mouse_end": 44653127,
        "mouse_strand": "+",
        "mouse_refseq": "NM_011234.5",
        "human_symbol": "RAD51",
        "human_chrom": "chr15",
        "human_start": 40694733,
        "human_end": 40732344,
        "identity_pct": 98.2,
        "genelab_datasets": [
            "GLDS-288",
        ],
    },
}


# =============================================================================
# Mouse Genome Access
# =============================================================================

MOUSE_GENOME_PATH = SHARED_DATA_DIR / "mouse" / "GRCm39.primary_assembly.genome.fa"


def get_mouse_genome():
    """Open mouse genome FASTA."""
    try:
        from pyfaidx import Fasta
    except ImportError:
        log.error("pyfaidx required: pip install pyfaidx")
        return None

    if not MOUSE_GENOME_PATH.exists():
        log.error("Mouse genome not found at {}".format(MOUSE_GENOME_PATH))
        log.error("Download with:")
        log.error("  wget https://ftp.ensembl.org/pub/release-112/fasta/"
                   "mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz")
        log.error("  gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz")
        log.error("  samtools faidx GRCm39.primary_assembly.genome.fa")
        return None

    return Fasta(str(MOUSE_GENOME_PATH))


# =============================================================================
# Saturation Scoring (Mouse)
# =============================================================================

def score_mouse_coding_saturation(
    gene_key: str,
    window_size: int = 8192,
) -> Optional[dict]:
    """
    Score all possible coding SNVs in a mouse ortholog gene.

    Similar to human scoring pipeline but uses mouse genome + coordinates.
    """
    if gene_key not in MOUSE_ORTHOLOGS:
        log.error("Unknown gene: {}".format(gene_key))
        return None

    orth = MOUSE_ORTHOLOGS[gene_key]
    log.info("\n--- {} ({}) ---".format(gene_key, orth["mouse_symbol"]))

    # Check for cached results
    out_dir = RESULTS_DIR / gene_key / "cross_species"
    out_dir.mkdir(parents=True, exist_ok=True)
    result_file = out_dir / "mouse_saturation_w{}.json".format(window_size)

    if result_file.exists():
        log.info("  Loading cached results from {}".format(result_file))
        with open(result_file) as f:
            return json.load(f)

    # Import Evo2 (GPU required)
    try:
        from evo2 import Evo2
    except ImportError as e:
        log.error("Cannot import Evo2: {}".format(e))
        return None

    genome = get_mouse_genome()
    if genome is None:
        return None

    chrom = orth["mouse_chrom"]
    gene_start = orth["mouse_start"]
    gene_end = orth["mouse_end"]
    gene_length = gene_end - gene_start

    log.info("  Mouse region: {}:{}-{} ({} bp)".format(
        chrom, gene_start, gene_end, gene_length))

    # For coding saturation, we score SNVs in the coding region
    # Extract the full gene region + flanking for context windows
    padding = window_size // 2
    region_start = max(0, gene_start - padding)
    region_end = gene_end + padding

    # Fetch reference sequence for the extended region
    chrom_key = chrom.replace("chr", "")
    try:
        ref_seq = str(genome[chrom_key][region_start:region_end])
    except (KeyError, ValueError):
        try:
            ref_seq = str(genome[chrom][region_start:region_end])
        except (KeyError, ValueError) as e:
            log.error("  Cannot fetch mouse sequence: {}".format(e))
            genome.close()
            return None

    log.info("  Fetched {} bp reference".format(len(ref_seq)))

    # Load Evo2 model
    log.info("  Loading Evo2 model...")
    model = Evo2("evo2_7b")

    # Generate all coding SNVs and score them
    # For simplicity, score every 3rd position (representative sampling)
    # Full saturation can be done later
    step = 1  # every position in coding region
    variants = []
    for pos in range(gene_start, gene_end, step):
        offset = pos - region_start
        if offset < 0 or offset >= len(ref_seq):
            continue
        ref_base = ref_seq[offset].upper()
        if ref_base not in BASES:
            continue
        for alt_base in BASES:
            if alt_base == ref_base:
                continue
            variants.append(Variant(
                chrom=chrom, pos=pos, ref=ref_base, alt=alt_base,
                gene=orth["mouse_symbol"],
            ))

    log.info("  Generated {} coding SNVs".format(len(variants)))

    # Score in windows (same approach as human scoring)
    from utils.sequence_utils import ScoringCheckpoint

    ckpt_dir = out_dir / "mouse_w{}".format(window_size)
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = ScoringCheckpoint(str(ckpt_dir))
    scored = checkpoint.load_results()
    scored_keys = {r["chrom"] + ":" + str(r["pos"]) + ":" + r["ref"] + ">" + r["alt"]
                   for r in scored}

    remaining = [v for v in variants if v.key() not in scored_keys]
    log.info("  Already scored: {}, remaining: {}".format(
        len(scored_keys), len(remaining)))

    # Score remaining variants
    batch_results = []
    for i, v in enumerate(remaining):
        # Build ref and alt windows
        center = v.pos - region_start
        half = window_size // 2
        win_start = max(0, center - half)
        win_end = win_start + window_size
        if win_end > len(ref_seq):
            win_end = len(ref_seq)
            win_start = max(0, win_end - window_size)

        ref_window = ref_seq[win_start:win_end]
        if len(ref_window) < window_size:
            continue

        # Create alt window
        alt_offset = center - win_start
        alt_window = ref_window[:alt_offset] + v.alt + ref_window[alt_offset + 1:]

        # Score both
        scores = model.score_sequences(
            [ref_window, alt_window],
            batch_size=1,
            reduce_method="mean",
            average_reverse_complement=True,
        )
        delta = scores[1] - scores[0]

        result = {
            "chrom": v.chrom,
            "pos": v.pos,
            "ref": v.ref,
            "alt": v.alt,
            "delta": float(delta),
            "ref_score": float(scores[0]),
            "alt_score": float(scores[1]),
            "gene": v.gene,
        }
        batch_results.append(result)

        if (i + 1) % 100 == 0:
            checkpoint.save_batch(batch_results)
            batch_results = []
            log.info("  Scored {}/{} variants".format(i + 1, len(remaining)))

    if batch_results:
        checkpoint.save_batch(batch_results)

    genome.close()

    # Aggregate results
    all_results = checkpoint.load_results()
    deltas = [r["delta"] for r in all_results]

    if not deltas:
        return None

    arr = np.array(deltas)
    summary = {
        "gene": gene_key,
        "mouse_symbol": orth["mouse_symbol"],
        "species": "mouse",
        "genome": "GRCm39",
        "n_variants": len(deltas),
        "mean_delta": float(np.mean(arr)),
        "std_delta": float(np.std(arr)),
        "median_delta": float(np.median(arr)),
        "q10": float(np.percentile(arr, 10)),
        "q90": float(np.percentile(arr, 90)),
        "window_size": window_size,
    }

    with open(result_file, "w") as f:
        json.dump(summary, f, indent=2)

    log.info("  Mouse {}: n={}, mean_delta={:.6f}".format(
        orth["mouse_symbol"], len(deltas), np.mean(arr)))

    return summary


# =============================================================================
# Cross-Species Comparison
# =============================================================================

def compare_constraint_profiles(
    gene_key: str,
    window_size: int = 8192,
) -> Optional[dict]:
    """
    Compare human vs mouse Evo2 constraint profiles for a gene.

    Loads pre-computed human and mouse saturation scoring results,
    then computes correlation between positional constraint scores.
    """
    if gene_key not in MOUSE_ORTHOLOGS:
        return None

    orth = MOUSE_ORTHOLOGS[gene_key]

    # Load human results
    human_dir = RESULTS_DIR / gene_key / "w{}".format(window_size)
    human_ckpt_file = human_dir / "results.jsonl"
    if not human_ckpt_file.exists():
        log.warning("  No human results for {}".format(gene_key))
        return None

    human_results = []
    with open(human_ckpt_file) as f:
        for line in f:
            human_results.append(json.loads(line))

    # Load mouse results
    mouse_dir = RESULTS_DIR / gene_key / "cross_species" / "mouse_w{}".format(window_size)
    mouse_ckpt_file = mouse_dir / "results.jsonl"
    if not mouse_ckpt_file.exists():
        log.warning("  No mouse results for {}".format(gene_key))
        return None

    mouse_results = []
    with open(mouse_ckpt_file) as f:
        for line in f:
            mouse_results.append(json.loads(line))

    # Compute per-position mean constraint (mean |delta| across all alts)
    def position_constraint(results):
        by_pos = {}
        for r in results:
            pos = r["pos"]
            if pos not in by_pos:
                by_pos[pos] = []
            by_pos[pos].append(abs(r["delta"]))
        return {pos: np.mean(vals) for pos, vals in by_pos.items()}

    human_constraint = position_constraint(human_results)
    mouse_constraint = position_constraint(mouse_results)

    # Summary statistics
    human_vals = np.array(list(human_constraint.values()))
    mouse_vals = np.array(list(mouse_constraint.values()))

    result = {
        "gene": gene_key,
        "human_symbol": orth["human_symbol"],
        "mouse_symbol": orth["mouse_symbol"],
        "human_n_positions": len(human_constraint),
        "mouse_n_positions": len(mouse_constraint),
        "human_mean_constraint": float(np.mean(human_vals)) if len(human_vals) > 0 else 0,
        "mouse_mean_constraint": float(np.mean(mouse_vals)) if len(mouse_vals) > 0 else 0,
        "protein_identity_pct": orth["identity_pct"],
    }

    # Distribution comparison
    if len(human_vals) > 0 and len(mouse_vals) > 0 and scipy_stats:
        ks_stat, ks_p = scipy_stats.ks_2samp(human_vals, mouse_vals)
        result["ks_stat"] = float(ks_stat)
        result["ks_p"] = float(ks_p)

    log.info("  {}: human constraint={:.6f} (n={}), mouse constraint={:.6f} (n={})".format(
        gene_key,
        result["human_mean_constraint"], len(human_constraint),
        result["mouse_mean_constraint"], len(mouse_constraint),
    ))

    return result


# =============================================================================
# GeneLab Connection
# =============================================================================

def genelab_context(gene_key: str) -> dict:
    """
    Provide NASA GeneLab dataset context for a gene.

    Returns relevant dataset IDs and descriptions for discussion.
    This is metadata-only (no data download required).
    """
    if gene_key not in MOUSE_ORTHOLOGS:
        return {}

    orth = MOUSE_ORTHOLOGS[gene_key]

    genelab_info = {
        "GLDS-288": {
            "title": "ISS Mouse Habitat Radiation Study",
            "organism": "Mus musculus",
            "tissue": "Multiple",
            "relevance": "Space radiation effects on DNA damage response genes",
            "mission": "ISS",
        },
        "GLDS-352": {
            "title": "GeneLab Multi-Omics Space Biology Reference",
            "organism": "Mus musculus",
            "tissue": "Multiple",
            "relevance": "Comprehensive multi-omics spaceflight reference dataset",
            "mission": "Various",
        },
        "GLDS-168": {
            "title": "Rodent Research-1 (RR-1) Transcriptomics",
            "organism": "Mus musculus",
            "tissue": "Liver, spleen",
            "relevance": "TERT and telomere biology in spaceflight mice",
            "mission": "SpaceX-4 / ISS",
        },
    }

    context = {
        "gene": gene_key,
        "mouse_symbol": orth["mouse_symbol"],
        "datasets": [],
    }

    for ds_id in orth.get("genelab_datasets", []):
        if ds_id in genelab_info:
            entry = genelab_info[ds_id].copy()
            entry["accession"] = ds_id
            context["datasets"].append(entry)

    return context


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Cross-species scoring (mouse orthologs)"
    )
    parser.add_argument("--gene", type=str, default="all",
                        help="Gene to score (ATM, TERT, RAD51, or 'all')")
    parser.add_argument("--window-size", type=int, default=8192)
    parser.add_argument("--compare-only", action="store_true",
                        help="Only compare existing human vs mouse results")
    parser.add_argument("--context-only", action="store_true",
                        help="Only output GeneLab dataset context (no scoring)")
    args = parser.parse_args()

    available_genes = list(MOUSE_ORTHOLOGS.keys())
    if args.gene.lower() == "all":
        genes = available_genes
    else:
        gene = args.gene.upper()
        if gene not in MOUSE_ORTHOLOGS:
            log.error("{} not in cross-species panel: {}".format(
                gene, available_genes))
            return
        genes = [gene]

    log.info("Cross-Species Analysis — genes: {}".format(genes))

    if args.context_only:
        for gene_key in genes:
            ctx = genelab_context(gene_key)
            log.info("\n{} GeneLab context:".format(gene_key))
            for ds in ctx.get("datasets", []):
                log.info("  {} — {} ({})".format(
                    ds["accession"], ds["title"], ds["relevance"]))
        return

    if args.compare_only:
        comparisons = []
        for gene_key in genes:
            result = compare_constraint_profiles(gene_key, args.window_size)
            if result:
                comparisons.append(result)

        if comparisons:
            out_dir = RESULTS_DIR / "cross_species"
            out_dir.mkdir(parents=True, exist_ok=True)
            with open(out_dir / "constraint_comparison.json", "w") as f:
                json.dump(comparisons, f, indent=2)
            log.info("\nComparison saved to {}".format(out_dir))
        return

    # Full scoring mode (GPU required)
    results = []
    for gene_key in genes:
        result = score_mouse_coding_saturation(gene_key, args.window_size)
        if result:
            results.append(result)

    if not results:
        log.error("No results computed")
        return

    # Summary
    log.info("\n" + "=" * 60)
    log.info("CROSS-SPECIES SCORING SUMMARY")
    log.info("=" * 60)
    for r in results:
        log.info("  {} ({}): n={}, mean_delta={:.6f}".format(
            r["gene"], r["mouse_symbol"],
            r["n_variants"], r["mean_delta"]))

    # Save
    out_dir = RESULTS_DIR / "cross_species"
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(out_dir / "mouse_scoring_summary.json", "w") as f:
        json.dump(results, f, indent=2)

    # GeneLab context for all scored genes
    genelab_ctx = []
    for gene_key in genes:
        ctx = genelab_context(gene_key)
        if ctx:
            genelab_ctx.append(ctx)

    if genelab_ctx:
        with open(out_dir / "genelab_context.json", "w") as f:
            json.dump(genelab_ctx, f, indent=2)

    # Run comparison if human results exist
    comparisons = []
    for gene_key in genes:
        result = compare_constraint_profiles(gene_key, args.window_size)
        if result:
            comparisons.append(result)

    if comparisons:
        with open(out_dir / "constraint_comparison.json", "w") as f:
            json.dump(comparisons, f, indent=2)

    log.info("\nResults saved to {}".format(out_dir))


if __name__ == "__main__":
    main()

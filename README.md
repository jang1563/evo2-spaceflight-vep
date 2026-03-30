# Evo2 Zero-Shot Variant Effect Prediction for Spaceflight Radiation-Response Genes

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Status**: Private repository. Will be made public upon manuscript acceptance.

## Overview

This repository contains the analysis pipeline for zero-shot variant effect prediction (VEP) across coding and non-coding regions of 10 spaceflight radiation-response genes using the [Evo2](https://github.com/ARC-Institute/evo2) genomic foundation model (7B parameters).

Evo2 scores sequence likelihood without any task-specific training, enabling variant effect assessment in regions inaccessible to protein-structure-based tools — particularly non-coding regulatory elements critical for genes like TERT, ATM, and NFE2L2.

**Key capabilities demonstrated:**
- Zero-shot VEP across coding AND non-coding regions from a single model
- Calibration against deep mutational scanning (DMS) data for 4 control genes
- Non-coding validation via TERT promoter MPRA and ENCODE regulatory element enrichment
- Functional impact profiling of radiation-characteristic mutations (C>A via 8-oxoG, microhomology deletions)
- PP3/BP4-level computational evidence for thousands of ClinVar VUS

## Gene Panel

| Gene | Type | Spaceflight Relevance |
|------|------|----------------------|
| **BRCA1** | DMS control | DNA repair; risk allele in Rutter et al. 2024 |
| **TP53** | DMS control | Most mutated gene in astronaut clonal hematopoiesis (7/14 astronauts) |
| **CHEK2** | DMS control | DNA damage checkpoint kinase |
| **DNMT3A** | DMS control | 2nd most mutated in astronaut CH; found in Inspiration4 crew |
| **TERT** | Novel target | Telomere biology; 14.5% elongation in NASA Twins Study |
| **ATM** | Novel target | DSB sensor for GCR/HZE radiation damage |
| **NFE2L2** | Novel target | Master antioxidant TF; 4 ISS mouse studies show protective role |
| **CLOCK** | Novel target | Circadian rhythm; 16 sunrises/day disruption on ISS |
| **MSTN** | Novel target | Muscle wasting protection in microgravity |
| **RAD51** | Novel target | Core HR recombinase for strand invasion at DSBs |

## Three-Layer Validation Framework

1. **DMS Calibration** — Spearman correlation with deep mutational scanning fitness scores (BRCA1, TP53, CHEK2, DNMT3A); gene-specific threshold calibration via Pejaver 2022 likelihood ratios
2. **ClinVar Validation** — Independent validation on ClinVar pathogenic/benign variants (>=2-star review), stratified by review status
3. **Non-Coding Validation** — TERT promoter MPRA (Kircher 2019), ENCODE regulatory element enrichment (permutation-based), comparison with LINSIGHT/Eigen/ncER/CADD

## Repository Structure

```
scripts/
  00_window_ablation.py       # Window size ablation study (4K-64K)
  01_prepare_variants.py      # Generate variant sequences (SNVs + indels)
  02_score_variants.py        # Core Evo2 scoring engine with checkpointing
  03_calibrate_dms.py         # DMS calibration + Pejaver LR thresholds
  04_validate_clinvar.py      # ClinVar validation + temporal split
  05_score_noncoding.py       # Non-coding scoring + TERT MPRA + ENCODE
  06_score_indels.py          # Indel scoring (in-frame vs frameshift)
  07_benchmark_tools.py       # Multi-tool comparison (AM, REVEL, CADD, etc.)
  08_radiation_signatures.py  # Radiation-characteristic mutation scoring
  09_cross_species.py         # Mouse ortholog scoring (GRCm39)
  10_entropy_landscape.py     # Positional entropy + ENCODE overlay
  11_astronaut_variants.py    # Published astronaut mutation scoring
  12_make_figures.py          # Publication figure generation
  benchmark_throughput.py     # GPU throughput benchmarking
  download_mavedb.py          # DMS data downloader (MaveDB API)
  utils/
    config.py                 # Centralized path configuration
    gene_coordinates.py       # Gene coordinates, exon boundaries (GRCh38)
    sequence_utils.py         # Window extraction, variant generation, ref dedup
    clinvar_parser.py         # ClinVar VCF parser with star-level stratification
    benchmarking.py           # AUROC/AUPRC/calibration/Pejaver LR computation
  slurm/
    submit_phase0.sh            # Phase 0: data download + setup
    submit_phase1_ablation.sh   # Phase 1: window ablation SLURM jobs
    template.sh                 # Generic SLURM job template
```

## Methods

### Scoring Formula

```
score(S) = (1/N) * SUM log P(s_{t+1} | s_1, ..., s_t)   # mean autoregressive log-prob
delta = score(S_alt) - score(S_ref)                        # per window
pathogenicity_priority = -delta                             # more negative = higher priority
```

All scores use reverse-complement averaging (`average_reverse_complement=True`) for strand-symmetric assessment.

### Window Size

Optimal window size selected via ablation across 4K, 8K, 16K, 32K, and 64K bp using the `evo2_7b` (<=8K) and `evo2_7b_262k` (>8K) model variants, evaluated by AUROC against ClinVar P/LP vs B/LB for control genes.

### Hardware

- Model: Evo2 7B (bfloat16, ~13.2 GB VRAM)
- GPUs: NVIDIA A40 (48 GB) and A100 (40/80 GB)
- Throughput: ~6.7 seconds/variant with RC averaging on A40
- Scoring is deterministic at bfloat16 precision

## Data Sources

| Dataset | Source | Build |
|---------|--------|-------|
| Reference genome | GRCh38 no-alt analysis set | hg38 |
| ClinVar | NCBI ClinVar VCF | hg38 |
| BRCA1 DMS | Findlay 2018 (MaveDB urn:mavedb:00000097-0-2) | — |
| TP53 DMS | Giacomelli 2018 (MaveDB urn:mavedb:00000068-b-1) | — |
| CHEK2 DMS | McCarthy-Leo 2024 (MaveDB urn:mavedb:00001203-a-1) | — |
| DNMT3A DMS | Garcia et al. 2025 (bioRxiv 10.1101/2025.09.24.678339) | — |
| TERT MPRA | Kircher 2019 (MaveDB urn:mavedb:00000031-b-1) | — |
| AlphaMissense | Cheng et al. 2023 | hg38 |
| CADD | v1.7 (Rentzsch et al. 2021) | hg38 |
| REVEL | v1.3 (Ioannidis et al. 2016) | hg38 |
| gnomAD constraint | v4.1 | hg38 |
| ENCODE | ChIP-seq peaks (DHS, H3K27ac, H3K4me1, CTCF) | hg38 |
| LINSIGHT / Eigen / ncER | Precomputed non-coding scores | hg19 (lifted over) |

## Installation

### Requirements

- Python >= 3.11
- CUDA-capable GPU with >= 48 GB VRAM (A40/A100)
- PyTorch >= 2.5 with CUDA 12.1+
- flash-attn >= 2.7

### Setup

```bash
conda create -n evo2 python=3.11 -y
conda activate evo2

# Core dependencies
pip install evo2 torch flash-attn

# Analysis dependencies
pip install biopython pandas numpy matplotlib seaborn scikit-learn \
    pysam pyfaidx pyranges openpyxl cyvcf2 gffutils jupyter
```

### Running

```bash
# 1. Download reference data
bash data/download_data.sh

# 2. Download DMS data from MaveDB
python scripts/download_mavedb.py

# 3. Run window ablation (requires GPU)
python scripts/00_window_ablation.py

# 4. Prepare variants
python scripts/01_prepare_variants.py

# 5. Score variants (SLURM recommended)
sbatch scripts/slurm/submit_phase1_ablation.sh
```

## Figures

| Figure | Description |
|--------|-------------|
| Fig 1 | Study overview: gene panel, hazard mapping, pipeline schematic |
| Fig 2 | Window size ablation: AUROC vs window size, score stability |
| Fig 3 | DMS calibration: Evo2 vs DMS fitness (4 control genes) |
| Fig 4 | ClinVar validation: ROC/PR curves, multi-tool AUROC comparison |
| Fig 5 | Per-gene constraint landscapes (all 10 genes) |
| Fig 6 | Non-coding validation: TERT MPRA, ENCODE enrichment, tool comparison |
| Fig 7 | Radiation functional impact: vulnerability maps, mutation spectra |
| Fig 8 | Spaceflight integration: astronaut variants, TERT promoter heatmap |

## Citation

If you use this code or data, please cite:

```bibtex
@article{kim2026evo2vep,
  title={Zero-shot variant effect prediction across coding and non-coding regions of spaceflight radiation-response genes using the Evo2 genomic foundation model},
  author={Kim, JangKeun and Mason, Christopher E.},
  journal={Genome Biology},
  year={2026}
}
```

And the Evo2 model:

```bibtex
@article{nguyen2026evo2,
  title={Sequence modeling and design from molecular to genome scale with Evo 2},
  author={Nguyen, Eric and others},
  journal={Nature},
  year={2026},
  doi={10.1038/s41586-026-10176-5}
}
```

## Related Work

- [Evo2](https://github.com/ARC-Institute/evo2) — Genomic foundation model (ARC Institute)
- [Rutter et al. 2024](https://doi.org/10.1038/s41467-024-50532-5) — Protective alleles and precision healthcare for spaceflight (*Nature Communications*)
- [SOMA Atlas](https://doi.org/10.1038/s41586-024-07639-y) — Space Omics and Medical Atlas (*Nature*)
- [Brojakowska et al. 2022](https://doi.org/10.1038/s42003-022-03777-z) — Astronaut clonal hematopoiesis (*Communications Biology*)

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

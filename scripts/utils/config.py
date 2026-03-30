"""
Project paths configuration.

Shared resources (genome, HF cache) live at the evo2/ level.
Project-specific data, results, logs live under project_spaceflight_vep/.

Set the EVO2_ROOT environment variable to your scratch directory, e.g.:
    export EVO2_ROOT=/scratch/$USER/evo2
"""

import os
from pathlib import Path

# Shared evo2 root (genome, HF models, tmp)
# Set EVO2_ROOT env var to your HPC scratch path, e.g. /scratch/$USER/evo2
EVO2_ROOT = Path(os.environ.get("EVO2_ROOT", "/path/to/your/scratch/evo2"))

# Project-specific directory
PROJECT_DIR = EVO2_ROOT / "project_spaceflight_vep"

# Shared data
SHARED_DATA_DIR = EVO2_ROOT / "data"
GENOME_PATH = SHARED_DATA_DIR / "GRCh38" / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Project-specific paths
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"
LOGS_DIR = PROJECT_DIR / "logs"
SCRIPTS_DIR = PROJECT_DIR / "scripts"

# ClinVar and annotation (shared — downloaded once, used by all projects)
CLINVAR_VCF = SHARED_DATA_DIR / "clinvar" / "clinvar.vcf.gz"
ANNOTATION_GFF = SHARED_DATA_DIR / "annotation" / "GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
ALPHAMISSENSE_TSV = SHARED_DATA_DIR / "benchmarks" / "alphamissense" / "AlphaMissense_hg38.tsv.gz"
GNOMAD_CONSTRAINT = SHARED_DATA_DIR / "gnomad" / "gnomad.v4.1.constraint_metrics.tsv"
LIFTOVER_CHAIN = SHARED_DATA_DIR / "liftover" / "hg19ToHg38.over.chain.gz"

# DMS, MPRA, ENCODE data (in shared area — single source of truth)
DMS_DIR = SHARED_DATA_DIR / "dms"
MPRA_DIR = SHARED_DATA_DIR / "mpra"
ENCODE_DIR = SHARED_DATA_DIR / "encode"

# HuggingFace cache
HF_HOME = EVO2_ROOT / "huggingface"

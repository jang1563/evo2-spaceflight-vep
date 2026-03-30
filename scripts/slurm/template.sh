#!/bin/bash
# =============================================================================
# SLURM Job Template for Evo2 VEP Scoring
# Cayuga HPC — scu-gpu partition (A40 48GB / A100 40-80GB)
#
# Usage (pass ALL resource parameters via sbatch command line):
#   sbatch --gres=gpu:1 --mem=64G --time=10:00:00 \
#     --export=SCRIPT=02_score_variants.py,GENE=BRCA1,WINDOW=8192 template.sh
#
# Env vars (passed via --export):
#   SCRIPT  - Python script to run (required)
#   GENE    - Gene symbol (optional, passed to script)
#   WINDOW  - Window size (optional, default 8192)
#
# NOTE: GPU count, memory, and time MUST be set via sbatch command-line flags,
# not env vars, because #SBATCH directives cannot expand shell variables.
# Set EVO2_ROOT before running, e.g.: export EVO2_ROOT=/scratch/$USER/evo2
# =============================================================================

#SBATCH --job-name=evo2_vep
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

echo "=== Evo2 VEP Job: $(date) ==="
echo "Host: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Script: ${SCRIPT:-unset}"
echo "Gene: ${GENE:-all}"
echo "Window: ${WINDOW:-8192}"
echo ""

# GPU info
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader

# Resolve EVO2_ROOT (set in environment or via --export)
EVO2_ROOT="${EVO2_ROOT:-/path/to/your/scratch/evo2}"

# Environment — update CONDA_BASE to your miniconda/conda path
CONDA_BASE="${CONDA_BASE:-/home/fs01/$USER/miniconda3/miniconda3}"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate evo2

export HF_HOME="${HF_HOME:-$EVO2_ROOT/../huggingface}"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TMPDIR="$EVO2_ROOT/tmp"
export PIP_CACHE_DIR="$EVO2_ROOT/tmp/pip_cache"

SCRIPTS_DIR="$EVO2_ROOT/project_spaceflight_vep/scripts"

# Run
echo ""
echo "=== Running ${SCRIPT} ==="
cd "$SCRIPTS_DIR"

python "${SCRIPT}" \
    ${GENE:+--gene "$GENE"} \
    ${WINDOW:+--window-size "$WINDOW"}

echo ""
echo "=== Done: $(date) ==="

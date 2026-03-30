#!/bin/bash
# =============================================================================
# Phase 1: Window size ablation study
# Submit SLURM jobs for each gene × window size combination
# Set EVO2_ROOT before running: export EVO2_ROOT=/scratch/$USER/evo2
# =============================================================================

set -euo pipefail

SBATCH="/opt/ohpc/pub/software/slurm/24.05.2/bin/sbatch"
EVO2_ROOT="${EVO2_ROOT:-/path/to/your/scratch/evo2}"
SCRIPTS_DIR="$EVO2_ROOT/project_spaceflight_vep/scripts"
LOGS_DIR="$EVO2_ROOT/project_spaceflight_vep/logs"
CONDA_BASE="${CONDA_BASE:-/home/fs01/$USER/miniconda3/miniconda3}"

echo "=== Phase 1: Window Ablation ==="

# Gene × window size matrix
GENES="BRCA1 TP53 CHEK2"

# Small windows (evo2_7b) — shorter time needed
for GENE in $GENES; do
    for WSIZE in 4096 8192; do
        JOB_NAME="ablation_${GENE}_w${WSIZE}"
        echo "Submitting $JOB_NAME..."
        $SBATCH \
            --job-name="$JOB_NAME" \
            --partition=scu-gpu \
            --gres=gpu:1 \
            --cpus-per-task=8 \
            --mem=64G \
            --time=04:00:00 \
            --output="$LOGS_DIR/${JOB_NAME}_%j.out" \
            --error="$LOGS_DIR/${JOB_NAME}_%j.err" \
            --wrap="
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate evo2
export HF_HOME=\${HF_HOME:-$(dirname $EVO2_ROOT)/huggingface}
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TMPDIR=$EVO2_ROOT/tmp
cd $SCRIPTS_DIR
python 00_window_ablation.py --gene $GENE --window-size $WSIZE
"
    done
done

# Larger windows (evo2_7b_262k) — more time, more memory
for GENE in $GENES; do
    for WSIZE in 16384 32768 65536; do
        JOB_NAME="ablation_${GENE}_w${WSIZE}"
        echo "Submitting $JOB_NAME..."
        $SBATCH \
            --job-name="$JOB_NAME" \
            --partition=scu-gpu \
            --gres=gpu:1 \
            --cpus-per-task=8 \
            --mem=96G \
            --time=08:00:00 \
            --output="$LOGS_DIR/${JOB_NAME}_%j.out" \
            --error="$LOGS_DIR/${JOB_NAME}_%j.err" \
            --wrap="
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate evo2
export HF_HOME=\${HF_HOME:-$(dirname $EVO2_ROOT)/huggingface}
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TMPDIR=$EVO2_ROOT/tmp
cd $SCRIPTS_DIR
python 00_window_ablation.py --gene $GENE --window-size $WSIZE
"
    done
done

echo ""
echo "Submitted 15 jobs (3 genes × 5 window sizes)"
echo "Monitor with: squeue -u \$USER"

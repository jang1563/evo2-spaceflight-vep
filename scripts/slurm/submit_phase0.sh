#!/bin/bash
# =============================================================================
# Phase 0: Infrastructure setup and data download
# Run on login node (no GPU needed)
# Set EVO2_ROOT before running: export EVO2_ROOT=/scratch/$USER/evo2
# =============================================================================

set -euo pipefail

SBATCH="/opt/ohpc/pub/software/slurm/24.05.2/bin/sbatch"
EVO2_ROOT="${EVO2_ROOT:-/path/to/your/scratch/evo2}"
SHARED_DATA_DIR="$EVO2_ROOT/data"
PROJECT_DIR="$EVO2_ROOT/project_spaceflight_vep"
PROJECT_DATA_DIR="$PROJECT_DIR/data"

echo "=== Phase 0: Setup ==="

# 1. Install missing packages
echo "Installing missing packages..."
CONDA_BASE="${CONDA_BASE:-/home/fs01/$USER/miniconda3/miniconda3}"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate evo2
export TMPDIR="$EVO2_ROOT/tmp"
export PIP_CACHE_DIR="$EVO2_ROOT/tmp/pip_cache"

pip install cyvcf2 gffutils requests 2>&1 | tail -3

# 2. Create shared data directories
echo "Creating shared data directories..."
mkdir -p "$SHARED_DATA_DIR"/{GRCh38,clinvar,annotation,benchmarks/{cadd,alphamissense,spliceai,revel},noncoding/{linsight,eigen,ncer},liftover,encode,gnomad,mouse}

# 3. Create project-specific data directories
echo "Creating project data directories..."
mkdir -p "$PROJECT_DATA_DIR"/{dms/{brca1,tp53,chek2,dnmt3a},mpra,encode}

# 4. Download data (run download script)
echo "Starting data downloads..."
bash "$EVO2_ROOT/data/download_data.sh"

# 5. Download evo2_7b_262k model (for window ablation)
echo ""
echo "=== Downloading evo2_7b_262k model ==="
export HF_HOME="${HF_HOME:-$EVO2_ROOT/../huggingface}"
python3 -c "
from huggingface_hub import snapshot_download
print('Downloading evo2_7b_262k...')
snapshot_download('arcinstitute/evo2_7b_262k')
print('Done.')
"

echo ""
echo "=== Phase 0 complete ==="

#!/bin/bash
#SBATCH --job-name=predraft
#SBATCH --output=/athena/masonlab/scratch/users/jak4013/evo2/project_spaceflight_vep/logs/predraft_%j.out
#SBATCH --error=/athena/masonlab/scratch/users/jak4013/evo2/project_spaceflight_vep/logs/predraft_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=scu-cpu

# Pre-draft analyses — no GPU needed (pure CPU analysis)
# Runs all 8 analyses: matched benchmark, DMS multitool, transversion control,
# ensemble, MPRA multitool, calibration regen, indel regen, variant count

set -euo pipefail

# Activate conda
source /home/fs01/jak4013/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate evo2

export EVO2_ROOT=/athena/masonlab/scratch/users/jak4013/evo2

cd /athena/masonlab/scratch/users/jak4013/evo2/project_spaceflight_vep/scripts

echo "Starting pre-draft analyses at $(date)"
echo "Python: $(which python)"

# Run all analyses
python 13_predraft_analyses.py --analysis all

echo "Completed at $(date)"

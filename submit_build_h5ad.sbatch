#!/bin/bash
#SBATCH --job-name=build_h5ad
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=build_h5ad_%j.out
#SBATCH --error=build_h5ad_%j.err

echo "Job started on $(hostname) at $(date)"

# Activate Conda environment
source /sc/kzd/home/abdullah/miniconda/etc/profile.d/conda.sh
conda activate bioinfo-notebook

# Run the conversion script
cd path/to/project_root
python3 02_build_h5ad.py

echo "Job completed at $(date)"

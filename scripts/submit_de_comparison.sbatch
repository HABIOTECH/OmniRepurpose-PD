#!/bin/bash
#SBATCH --job-name=de_compare
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=de_compare_%j.out
#SBATCH --error=de_compare_%j.err

echo "=== Pipeline started at $(date) ==="
source /sc/kzd/home/abdullah/miniconda/etc/profile.d/conda.sh

cd /sc/kzd/home/abdullah/bioinformatics/PK-multiomics

echo "--- [1/3] Running R DESeq2 (deseq2_env) ---"
conda activate deseq2_env
Rscript scripts/02c_run_deseq2.R

echo "--- [2/3] Running PyDESeq2 (base) ---"
conda activate base
python3 scripts/02d_run_pydeseq2.py

echo "--- [3/3] Running Comparison Script (base) ---"
python3 scripts/02e_compare_de.py

echo "=== Pipeline finished at $(date) ==="

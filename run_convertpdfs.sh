#!/bin/bash
#SBATCH --array=1-27      # adjust to number of biomarkers
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/convert_%A_%a.out

module load python/3.11   # or your cluster’s Python module

BIOMARKER=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /work/users/m/c/mcgeet/DMETpaper/DMET-results/biomarkernames_pyth.txt | tr -d '"')

INPUT_DIR="/work/users/m/c/mcgeet/DMETpaper/figs/${BIOMARKER}/rint"
OUTPUT_DIR="/work/users/m/c/mcgeet/DMETpaper/figures/${BIOMARKER}/jpgs/rint"

python /work/users/m/c/mcgeet/DMETpaper/DMET-results/convert_pdfs.py "$INPUT_DIR" "$OUTPUT_DIR"

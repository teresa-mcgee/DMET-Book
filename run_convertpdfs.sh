#!/bin/bash
#SBATCH --array=1-27      # adjust to number of biomarkers
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/convert_%A_%a.out

module load python/3.11   # or your cluster’s Python module

BIOMARKER=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /work/users/m/c/mcgeet/DMETpaper/DMET-results/biomarkernames_pyth.txt | tr -d '"')

# NOTE: this predates the current figure pipeline (run_final_fig.sh -> R/dev/final_figure_generation.R),
# which writes finished jpgs directly into DMET-results/figures/. Kept for reference; INPUT_DIR below
# now points at the archived source tree since figs/ was moved to _archive/figs/ in the 2026 cleanup.
INPUT_DIR="/work/users/m/c/mcgeet/DMETpaper/_archive/figs/${BIOMARKER}/rint"
OUTPUT_DIR="/work/users/m/c/mcgeet/DMETpaper/DMET-results/figures/${BIOMARKER}/jpgs/rint"

python /work/users/m/c/mcgeet/DMETpaper/DMET-results/convert_pdfs.py "$INPUT_DIR" "$OUTPUT_DIR"

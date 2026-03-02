#!/usr/bin/env bash
# slurm_bslmm_summarize.sh
# Runs summarize_bslmm.py after all BSLMM array jobs complete.
#
# Submit with a dependency on the array job:
#   sbatch --dependency=afterok:<ARRAY_JOBID> slurm_bslmm_summarize.sh
#
# Or submit manually once the array jobs are finished:
#   sbatch slurm_bslmm_summarize.sh

#SBATCH --job-name=bslmm_summarize
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --partition=general            # ← change to your cluster's partition name
#SBATCH --output=logs/bslmm_summarize_%j.out
#SBATCH --error=logs/bslmm_summarize_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your@email.edu     # ← replace with your address

set -euo pipefail
mkdir -p logs

# Load Python environment — adjust for your cluster
# e.g. module load python/3.11, conda activate myenv, or source venv/bin/activate
module load python 2>/dev/null || true

BASE="$(dirname "$0")"

echo "[$(date '+%H:%M:%S')] Running BSLMM summary..."
python3 "${BASE}/summarize_bslmm.py"
echo "[$(date '+%H:%M:%S')] Summary complete. Results in ${BASE}/gemma_output/bslmm/"

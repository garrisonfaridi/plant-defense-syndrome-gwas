#!/usr/bin/env bash
# slurm_bslmm_array.sh
# SLURM array job: runs BSLMM for PC1–PC8 in parallel (one job per PC).
#
# Submit from the repo root:
#   sbatch slurm_bslmm_array.sh
#
# Once complete, submit the summary step:
#   sbatch --dependency=afterok:<JOBID> slurm_bslmm_summarize.sh
#
# ── Adjust the #SBATCH directives below for your cluster ──────────────────

#SBATCH --job-name=bslmm_gwas
#SBATCH --array=1-8                    # one task per PC
#SBATCH --cpus-per-task=8              # GEMMA uses OpenMP threads
#SBATCH --mem=24G                      # ~2M SNPs; 16G usually sufficient, 24G safe
#SBATCH --time=24:00:00                # 10M samples ~ 13–16 h; 24 h gives headroom
#SBATCH --partition=general            # ← change to your cluster's partition name
#SBATCH --output=logs/bslmm_PC%a_%j.out
#SBATCH --error=logs/bslmm_PC%a_%j.err
#SBATCH --mail-type=END,FAIL           # ← comment out if not wanted
#SBATCH --mail-user=your@email.edu     # ← replace with your address

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
set -euo pipefail
mkdir -p logs

# Load GEMMA — adjust module name for your HPC environment
# Common alternatives: module load gemma, module load GEMMA/0.98.5
module load gemma 2>/dev/null || true

# If GEMMA is not in a module, set an explicit path here:
# GEMMA_BIN=/path/to/gemma
GEMMA_BIN=${GEMMA_BIN:-gemma}

# Tell GEMMA how many threads to use (matches --cpus-per-task)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

# ---------------------------------------------------------------------------
# Paths  ── all relative to the repo root; adjust BASE if needed
# ---------------------------------------------------------------------------
BASE="$(dirname "$0")"                 # directory containing this script
GENO="${BASE}/genotypes/genotypes_305"
PHENO="${BASE}/phenotype.txt"
KINSHIP="${BASE}/gemma_output/kinship/kinship.cXX.txt"
BSLMM_DIR="${BASE}/gemma_output/bslmm"

mkdir -p "${BSLMM_DIR}"

# Guard: kinship matrix must exist (compute with run_lmm.sh or run_gemma.sh first)
if [[ ! -f "${KINSHIP}" ]]; then
  echo "ERROR: kinship matrix not found at ${KINSHIP}" >&2
  echo "       Run run_lmm.sh first to generate it." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# BSLMM — publication-quality parameters
#   -w  500000   burn-in steps
#   -s 10000000  sampling steps  (10× exploratory run)
#   -rpace 10    record every 10th sample (reduces .hyp/.gamma file sizes)
# ---------------------------------------------------------------------------
PC=${SLURM_ARRAY_TASK_ID}

echo "[$(date '+%H:%M:%S')] Starting BSLMM PC${PC} on $(hostname)"
echo "  OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "  burn-in=500000  samples=10000000"

${GEMMA_BIN} \
  -bfile  "${GENO}"     \
  -p      "${PHENO}"    \
  -k      "${KINSHIP}"  \
  -n      "${PC}"       \
  -bslmm  1             \
  -w      500000        \
  -s      10000000      \
  -rpace  10            \
  -o      "bslmm_PC${PC}" \
  -outdir "${BSLMM_DIR}"

echo "[$(date '+%H:%M:%S')] BSLMM PC${PC} complete"

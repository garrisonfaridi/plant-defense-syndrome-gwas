#!/usr/bin/env bash
# run_lmm.sh — Kinship matrix + LMM GWAS + Manhattan/QQ plots for 8 pseudotraits
# Usage: bash scripts/shell/run_lmm.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
GENO="${REPO_ROOT}/genotypes/genotypes_305"
PHENO="${REPO_ROOT}/data/phenotype.txt"
KINSHIP="${REPO_ROOT}/results/kinship/kinship.cXX.txt"
LMM_DIR="${REPO_ROOT}/results/lmm"

ts() { echo "[$(date '+%H:%M:%S')] $*"; }

# ---------------------------------------------------------------------------
# Step 1 — Compute centered relatedness (kinship) matrix
# ---------------------------------------------------------------------------
mkdir -p "${REPO_ROOT}/results/kinship"
ts "Computing kinship matrix..."
gemma \
  -bfile "${GENO}" \
  -p "${PHENO}" \
  -gk 1 \
  -o kinship \
  -outdir "${REPO_ROOT}/results/kinship"
ts "Kinship done -> ${KINSHIP}"

# ---------------------------------------------------------------------------
# Step 2 — LMM for traits 1–8
# ---------------------------------------------------------------------------
mkdir -p "${LMM_DIR}"
COMMON_ARGS="-bfile ${GENO} -p ${PHENO} -k ${KINSHIP}"
ts "Starting LMM runs..."
for i in {1..8}; do
  ts "  LMM PC${i}..."
  gemma ${COMMON_ARGS} -n "${i}" -lmm 4 -o "lmm_PC${i}" -outdir "${LMM_DIR}"
  ts "  LMM PC${i} done"
done
ts "All LMM runs complete."

# ---------------------------------------------------------------------------
# Step 3 — Manhattan + QQ plots and summary statistics
# ---------------------------------------------------------------------------
ts "Generating Manhattan and QQ plots..."
python3 "${REPO_ROOT}/scripts/python/plot_lmm.py"
ts "Plotting done. Results in ${LMM_DIR}/"

#!/usr/bin/env bash
# run_lmm.sh — Kinship matrix + LMM GWAS + Manhattan plots for 8 pseudotraits
# Usage: bash Cirrone_Lab/run_lmm.sh
# Run from: /Users/garrisonfaridi/Documents/Documents_Cloud/BioCode/

set -euo pipefail

BASE="Cirrone_Lab"
GENO="${BASE}/genotypes/genotypes_305"
PHENO="${BASE}/phenotype.txt"
KINSHIP="${BASE}/gemma_output/kinship/kinship.cXX.txt"
LMM_DIR="${BASE}/gemma_output/lmm"

ts() { echo "[$(date '+%H:%M:%S')] $*"; }

# ---------------------------------------------------------------------------
# Step 1 — Compute centered relatedness (kinship) matrix
# ---------------------------------------------------------------------------
mkdir -p "${BASE}/gemma_output/kinship"
ts "Computing kinship matrix..."
gemma \
  -bfile "${GENO}" \
  -p "${PHENO}" \
  -gk 1 \
  -o kinship \
  -outdir "${BASE}/gemma_output/kinship"
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
# Step 3 — Manhattan plots and summary statistics
# ---------------------------------------------------------------------------
ts "Generating Manhattan plots..."
python3 "${BASE}/plot_lmm.py"
ts "Plotting done. Results in ${LMM_DIR}/"

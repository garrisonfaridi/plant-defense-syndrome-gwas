#!/usr/bin/env bash
# run_bslmm.sh — BSLMM GWAS + posterior summaries for 8 pseudotraits
# Requires kinship.cXX.txt to already exist (run run_lmm.sh first).
# Usage: bash scripts/shell/run_bslmm.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
GENO="${REPO_ROOT}/genotypes/genotypes_305"
PHENO="${REPO_ROOT}/data/phenotype.txt"
KINSHIP="${REPO_ROOT}/results/kinship/kinship.cXX.txt"
BSLMM_DIR="${REPO_ROOT}/results/bslmm"

ts() { echo "[$(date '+%H:%M:%S')] $*"; }

# Guard: kinship must exist before BSLMM can run.
if [[ ! -f "${KINSHIP}" ]]; then
  echo "ERROR: ${KINSHIP} not found. Run run_lmm.sh first to compute the kinship matrix."
  exit 1
fi

# ---------------------------------------------------------------------------
# Step 1 — BSLMM for traits 1–8  (100K burn-in, 1M samples)
# ---------------------------------------------------------------------------
mkdir -p "${BSLMM_DIR}"
COMMON_ARGS="-bfile ${GENO} -p ${PHENO} -k ${KINSHIP}"
ts "Starting BSLMM runs..."
for i in {1..8}; do
  ts "  BSLMM PC${i}..."
  gemma ${COMMON_ARGS} -n "${i}" -bslmm 1 -w 100000 -s 1000000 -o "bslmm_PC${i}" -outdir "${BSLMM_DIR}"
  ts "  BSLMM PC${i} done"
done
ts "All BSLMM runs complete."

# ---------------------------------------------------------------------------
# Step 2 — Posterior summaries and top SNP tables
# ---------------------------------------------------------------------------
ts "Summarizing BSLMM results..."
python3 "${REPO_ROOT}/scripts/python/summarize_bslmm.py"
ts "Done. Results in ${BSLMM_DIR}/"

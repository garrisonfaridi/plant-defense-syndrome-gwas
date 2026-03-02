#!/usr/bin/env bash
# run_gemma.sh — Full pipeline: kinship + LMM + BSLMM on 8 pseudotraits
# Delegates to run_lmm.sh then run_bslmm.sh sequentially.
# Usage: bash scripts/shell/run_gemma.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
bash "${SCRIPT_DIR}/run_lmm.sh"
bash "${SCRIPT_DIR}/run_bslmm.sh"

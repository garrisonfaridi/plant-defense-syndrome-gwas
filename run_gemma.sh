#!/usr/bin/env bash
# run_gemma.sh — Full pipeline: kinship + LMM + BSLMM on 8 pseudotraits
# Delegates to run_lmm.sh then run_bslmm.sh sequentially.
# Usage: bash Cirrone_Lab/run_gemma.sh
# Run from: /Users/garrisonfaridi/Documents/Documents_Cloud/BioCode/

set -euo pipefail

SCRIPT_DIR="$(dirname "$0")"
bash "${SCRIPT_DIR}/run_lmm.sh"
bash "${SCRIPT_DIR}/run_bslmm.sh"

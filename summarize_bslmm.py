"""
summarize_bslmm.py — BSLMM posterior summaries for 8 pseudotraits.

Run from: /Users/garrisonfaridi/Documents/Documents_Cloud/BioCode/
Outputs:
  Cirrone_Lab/gemma_output/bslmm/summary_hyp.csv
  Cirrone_Lab/gemma_output/bslmm/summary_bslmm.csv
  Cirrone_Lab/gemma_output/bslmm/top_snps_PC{i}.csv   (i = 1..8)
  Cirrone_Lab/gemma_output/bslmm/interpretation_bslmm.txt
"""

import os
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = "/Users/garrisonfaridi/Documents/Documents_Cloud/BioCode/Cirrone_Lab"
BSLMM_DIR = os.path.join(BASE, "gemma_output", "bslmm")
os.makedirs(BSLMM_DIR, exist_ok=True)

HYP_COLS = ["pve", "pge", "rho", "pi", "n_gamma"]

# ---------------------------------------------------------------------------
# Accumulators
# ---------------------------------------------------------------------------
hyp_rows = []       # one row per trait × hyperparameter
summ_rows = []      # one row per trait (credible-set summary)
interp_lines = []   # text interpretation

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
for i in range(1, 9):
    hyp_file = os.path.join(BSLMM_DIR, f"bslmm_PC{i}.hyp.txt")
    param_file = os.path.join(BSLMM_DIR, f"bslmm_PC{i}.param.txt")

    # Initialize fallback values; overwritten below if files exist.
    pve_med = pve_lo = pve_hi = np.nan
    pge_med = pge_lo = pge_hi = np.nan
    n_gamma_med = np.nan
    n_pip001 = np.nan

    # -----------------------------------------------------------------------
    # .hyp.txt — posterior samples of hyperparameters
    # -----------------------------------------------------------------------
    if not os.path.exists(hyp_file):
        print(f"[WARN] {hyp_file} not found — skipping PC{i} hyp")
    else:
        hyp = pd.read_csv(hyp_file, sep="\t", index_col=False)
        hyp.columns = hyp.columns.str.strip()
        hyp = hyp.dropna(axis=1, how="all")  # drop trailing all-NaN column from trailing tab
        # Align to expected columns (GEMMA may include extra cols)
        for col in HYP_COLS:
            if col not in hyp.columns:
                hyp[col] = np.nan

        # Single pass: compute stats and cache them to avoid re-querying the DataFrame.
        stats_cache = {}
        for col in HYP_COLS:
            vals = hyp[col].dropna().values
            med = float(np.median(vals))
            lo = float(np.percentile(vals, 2.5))
            hi = float(np.percentile(vals, 97.5))
            stats_cache[col] = (med, lo, hi)
            hyp_rows.append({
                "Trait": f"PC{i}",
                "Parameter": col,
                "Median": med,
                "CI_2.5": lo,
                "CI_97.5": hi,
            })

        pve_med, pve_lo, pve_hi = stats_cache["pve"]
        pge_med, pge_lo, pge_hi = stats_cache["pge"]
        n_gamma_med = stats_cache["n_gamma"][0]
        print(f"[PC{i}] PVE={pve_med:.3f} [{pve_lo:.3f},{pve_hi:.3f}]  "
              f"PGE={pge_med:.3f} [{pge_lo:.3f},{pge_hi:.3f}]  "
              f"n_gamma(median)={n_gamma_med:.1f}")

    # -----------------------------------------------------------------------
    # .param.txt — per-SNP posterior estimates
    # -----------------------------------------------------------------------
    if not os.path.exists(param_file):
        print(f"[WARN] {param_file} not found — skipping PC{i} param")
    else:
        param = pd.read_csv(param_file, sep="\t")
        # PIP = gamma column (posterior mean of inclusion indicator)
        param["pip"] = param["gamma"].astype(float)
        param["eff"] = param["alpha"].astype(float) * param["pip"]

        n_pip001 = int((param["pip"] > 0.01).sum())

        # Top 20 SNPs by PIP — write directly without retaining the intermediate variable.
        top_csv = os.path.join(BSLMM_DIR, f"top_snps_PC{i}.csv")
        top_snps = param.nlargest(20, "pip")[["chr", "rs", "ps", "n_miss", "alpha", "beta", "pip", "eff"]].reset_index(drop=True)
        top_snps.insert(0, "Trait", f"PC{i}")
        top_snps.to_csv(top_csv, index=False)
        print(f"  Top SNPs saved -> {top_csv}  (n_PIP>0.01: {n_pip001})")

    # -----------------------------------------------------------------------
    # Summary row
    # -----------------------------------------------------------------------
    summ_rows.append({
        "Trait": f"PC{i}",
        "PVE_median": pve_med,
        "PVE_2.5": pve_lo,
        "PVE_97.5": pve_hi,
        "PGE_median": pge_med,
        "PGE_2.5": pge_lo,
        "PGE_97.5": pge_hi,
        "n_gamma_median": n_gamma_med,
        "n_snps_PIP_gt_0.01": n_pip001,
    })

    # -----------------------------------------------------------------------
    # Interpretation
    # -----------------------------------------------------------------------
    if np.isnan(pve_med):
        interp_lines.append(f"PC{i}: Results not available (output files missing).")
        continue

    ci_excludes_zero = pve_lo > 0
    signal_str = (
        "95% CI for PVE excludes 0 — strong genetic signal detected."
        if ci_excludes_zero
        else "95% CI for PVE overlaps 0 — weak or no detectable genetic signal."
    )

    if not np.isnan(pge_med):
        arch_str = (
            f"PGE={pge_med:.3f} > 0.5 → most genetic variance comes from a small number "
            f"of large-effect loci (oligogenic architecture)."
            if pge_med > 0.5
            else f"PGE={pge_med:.3f} < 0.5 → genetic variance is spread across many small-effect "
                 f"SNPs (highly polygenic architecture)."
        )
    else:
        arch_str = "PGE not available."

    # Use filter(None, ...) to drop empty fragments and avoid double-spaces.
    ngam_str = (
        f"Estimated n_gamma (median)={n_gamma_med:.1f} causal SNPs contributing to genetic variance."
        if not np.isnan(n_gamma_med) else ""
    )
    pip_str = (
        f"{n_pip001} SNPs have PIP > 0.01 (candidates for follow-up)."
        if not np.isnan(n_pip001) else ""
    )

    interp_lines.append(" ".join(filter(None, [
        f"PC{i}: PVE={pve_med:.3f} [{pve_lo:.3f}, {pve_hi:.3f}] — "
        f"proportion of phenotypic variance explained (narrow-sense heritability estimate).",
        signal_str,
        arch_str,
        ngam_str,
        pip_str,
    ])))

# ---------------------------------------------------------------------------
# Write summary_hyp.csv
# ---------------------------------------------------------------------------
if hyp_rows:
    pd.DataFrame(hyp_rows).to_csv(os.path.join(BSLMM_DIR, "summary_hyp.csv"), index=False)
    print(f"\nHyperparameter summary -> {os.path.join(BSLMM_DIR, 'summary_hyp.csv')}")

# ---------------------------------------------------------------------------
# Write summary_bslmm.csv
# ---------------------------------------------------------------------------
if summ_rows:
    summ_df = pd.DataFrame(summ_rows)
    summ_csv = os.path.join(BSLMM_DIR, "summary_bslmm.csv")
    summ_df.to_csv(summ_csv, index=False)
    print(f"Credible-set summary -> {summ_csv}")
    print(summ_df.to_string(index=False))

# ---------------------------------------------------------------------------
# Write interpretation_bslmm.txt
# ---------------------------------------------------------------------------
interp_txt = os.path.join(BSLMM_DIR, "interpretation_bslmm.txt")
with open(interp_txt, "w") as fh:
    fh.write("BSLMM Interpretation — PC1 to PC8\n")
    fh.write("=" * 60 + "\n\n")
    fh.write(
        "Key parameters:\n"
        "  PVE   : Proportion of phenotypic variance explained by all SNPs (narrow-sense h²).\n"
        "  PGE   : Proportion of genetic variance from sparse (large-effect) SNPs vs polygenic background.\n"
        "  n_gamma: Effective number of SNPs with non-zero sparse effects.\n"
        "  PIP   : Posterior inclusion probability per SNP (>0.1 = candidate; >0.5 = strong).\n\n"
    )
    for line in interp_lines:
        fh.write(line + "\n\n")
print(f"Interpretation -> {interp_txt}")

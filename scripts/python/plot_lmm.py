"""
plot_lmm.py — Manhattan plots, QQ plots, and summary statistics for 8 LMM GWAS results.

Run from any directory:
    python scripts/python/plot_lmm.py

Outputs (relative to repo root):
  results/lmm/manhattan_PC{i}.png   (i = 1..8)
  results/lmm/qq_PC{i}.png          (i = 1..8)
  results/lmm/qq_panel.png          (2×4 merged panel)
  results/lmm/summary_lmm.csv
  results/lmm/interpretation_lmm.txt
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# ---------------------------------------------------------------------------
# Paths  (scripts/python/ → scripts/ → repo root)
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LMM_DIR = os.path.join(REPO_ROOT, "results", "lmm")
os.makedirs(LMM_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Chromosome layout constants (TAIR10, Arabidopsis thaliana)
# Computed once — reused across all 8 trait plots.
# ---------------------------------------------------------------------------
CHROM_SIZES = {1: 30427671, 2: 19698289, 3: 23459830, 4: 18585056, 5: 26975502}
CHROMS = [1, 2, 3, 4, 5]
CHR_COLORS = ["steelblue", "lightsteelblue"]

# Build chromosome offsets and midpoints once
_running = 0
CHR_OFFSETS = {}
for _c in CHROMS:
    CHR_OFFSETS[_c] = _running
    _running += CHROM_SIZES[_c]
GENOME_LEN = _running
CHR_MIDS = {c: CHR_OFFSETS[c] + CHROM_SIZES[c] / 2 for c in CHROMS}

# ---------------------------------------------------------------------------
# QQ plot helper
# ---------------------------------------------------------------------------
def _qq_arrays(pvals: np.ndarray, bonferroni: float):
    """Compute all arrays needed for a QQ plot. Returns a dict of precomputed values."""
    pvals = pvals[np.isfinite(pvals) & (pvals > 0)]
    pvals = np.sort(pvals)
    n = len(pvals)
    expected = np.arange(1, n + 1) / (n + 1)
    obs_log  = -np.log10(pvals)
    exp_log  = -np.log10(expected)
    chi2_obs  = stats.chi2.ppf(1 - pvals, df=1)
    lambda_gc = np.median(chi2_obs) / stats.chi2.ppf(0.5, df=1)
    ci_lo = -np.log10(stats.beta.ppf(0.975, np.arange(1, n + 1), np.arange(n, 0, -1)))
    ci_hi = -np.log10(stats.beta.ppf(0.025, np.arange(1, n + 1), np.arange(n, 0, -1)))
    sig_mask = pvals < bonferroni
    return dict(obs_log=obs_log, exp_log=exp_log, ci_lo=ci_lo, ci_hi=ci_hi,
                sig_mask=sig_mask, lambda_gc=lambda_gc)


def _draw_qq_ax(ax, arrays: dict, title: str, fontsize: int = 9) -> None:
    """Draw a QQ plot into an existing Axes object."""
    obs_log  = arrays["obs_log"]
    exp_log  = arrays["exp_log"]
    ci_lo    = arrays["ci_lo"]
    ci_hi    = arrays["ci_hi"]
    sig_mask = arrays["sig_mask"]
    lambda_gc = arrays["lambda_gc"]

    ax.fill_between(exp_log, ci_lo, ci_hi, color="lightsteelblue", alpha=0.5)
    ax.scatter(exp_log[~sig_mask], obs_log[~sig_mask],
               c="steelblue", s=3, alpha=0.6, linewidths=0, rasterized=True)
    if sig_mask.any():
        ax.scatter(exp_log[sig_mask], obs_log[sig_mask],
                   c="orangered", s=8, alpha=0.9, linewidths=0, zorder=5)
    lim = max(exp_log.max(), obs_log.max()) * 1.05
    ax.plot([0, lim], [0, lim], color="crimson", linewidth=0.9, linestyle="--")
    ax.set_xlim(0, exp_log.max() * 1.05)
    ax.set_ylim(0, obs_log.max() * 1.05)
    ax.set_title(f"{title}\n$\\lambda_{{GC}}={lambda_gc:.3f}$", fontsize=fontsize)
    ax.set_xlabel(r"Expected $-\log_{10}(p)$", fontsize=fontsize - 1)
    ax.set_ylabel(r"Observed $-\log_{10}(p)$", fontsize=fontsize - 1)
    ax.tick_params(labelsize=fontsize - 2)


def plot_qq(pvals: np.ndarray, title: str, outpath: str, bonferroni: float) -> float:
    """
    Save a standalone QQ plot and return λ_GC.
    Confidence band from Beta order-statistic distribution.
    """
    arrays = _qq_arrays(pvals, bonferroni)
    fig, ax = plt.subplots(figsize=(6, 6))
    _draw_qq_ax(ax, arrays, title, fontsize=11)
    # Add legend only on standalone plot
    import matplotlib.patches as mpatches
    ax.legend(
        handles=[
            mpatches.Patch(color="lightsteelblue", alpha=0.5, label="95% CI"),
            plt.Line2D([0], [0], color="crimson", linestyle="--", label="Expected (null)"),
            plt.scatter([], [], c="orangered", s=10, label="Significant SNPs"),
        ],
        fontsize=8, loc="upper left",
    )
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return arrays["lambda_gc"]


# ---------------------------------------------------------------------------
# Summary accumulators
# ---------------------------------------------------------------------------
summary_rows = []
interp_lines = []
qq_data = {}   # {pc_index: arrays dict} — collected for the merged panel

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
for i in range(1, 9):
    assoc_file = os.path.join(LMM_DIR, f"lmm_PC{i}.assoc.txt")
    if not os.path.exists(assoc_file):
        print(f"[WARN] {assoc_file} not found — skipping PC{i}")
        continue

    df = pd.read_csv(assoc_file, sep="\t")

    # Filter to autosomes 1–5 and drop rows with missing p-values
    df = df[df["chr"].isin(CHROMS)].copy()
    df = df.dropna(subset=["p_wald"])
    df["p_wald"] = df["p_wald"].astype(float)
    df["beta"] = df["beta"].astype(float)

    n_snps = len(df)
    bonferroni = 0.05 / n_snps
    bonf_log = -np.log10(bonferroni)
    df["neg_log10_p"] = -np.log10(df["p_wald"].clip(lower=1e-300))

    # Vectorized cumulative position — map offsets dict onto chr column
    df["cum_pos"] = df["chr"].astype(int).map(CHR_OFFSETS).fillna(0) + df["ps"].astype(int)

    # Significance mask computed once per trait
    sig_mask = df["p_wald"] < bonferroni
    sig = df[sig_mask]
    n_sig = len(sig)

    # Top SNP
    top = sig.loc[sig["p_wald"].idxmin()] if n_sig > 0 else df.loc[df["p_wald"].idxmin()]
    top_rs = top.get("rs", "NA")
    top_chr = int(top["chr"])
    top_pos = int(top["ps"])
    top_pval = top["p_wald"]
    top_beta = top["beta"]

    # Pre-split DataFrame by chromosome to avoid per-chromosome full scans
    chr_groups = dict(list(df.groupby("chr")))

    # -------------------------------------------------------------------
    # Manhattan plot
    # -------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(16, 5))
    sig_label_added = False

    for ci, chrom in enumerate(CHROMS):
        sub = chr_groups.get(chrom)
        if sub is None or sub.empty:
            continue
        color = CHR_COLORS[ci % 2]
        chrom_sig_mask = sub["p_wald"] < bonferroni
        non_sig_sub = sub[~chrom_sig_mask]
        sig_sub = sub[chrom_sig_mask]

        ax.scatter(
            non_sig_sub["cum_pos"],
            non_sig_sub["neg_log10_p"],
            c=color, s=4, alpha=0.7, linewidths=0, rasterized=True,
        )
        if not sig_sub.empty:
            label = "Significant SNPs" if not sig_label_added else ""
            ax.scatter(
                sig_sub["cum_pos"],
                sig_sub["neg_log10_p"],
                c="orangered", s=14, alpha=0.95, linewidths=0, zorder=5, label=label,
            )
            sig_label_added = True

    ax.axhline(bonf_log, color="crimson", linestyle="--", linewidth=1.2,
               label=f"Bonferroni threshold (p={bonferroni:.2e})")

    ax.set_xticks([CHR_MIDS[c] for c in CHROMS])
    ax.set_xticklabels([f"Chr{c}" for c in CHROMS], fontsize=10)
    ax.set_xlim(0, GENOME_LEN)
    ax.set_xlabel("Chromosome", fontsize=11)
    ax.set_ylabel(r"$-\log_{10}(p_{\rm wald})$", fontsize=11)
    ax.set_title(f"LMM GWAS — PC{i}  (n_sig={n_sig}, Bonf. p<{bonferroni:.2e})", fontsize=13)
    ax.legend(fontsize=9, loc="upper right")

    plt.tight_layout()
    out_png = os.path.join(LMM_DIR, f"manhattan_PC{i}.png")
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"[PC{i}] Manhattan saved -> {out_png}")

    # -------------------------------------------------------------------
    # QQ plot
    # -------------------------------------------------------------------
    arrays = _qq_arrays(df["p_wald"].values, bonferroni)
    qq_data[i] = arrays          # cache for merged panel
    lambda_gc = arrays["lambda_gc"]

    qq_png = os.path.join(LMM_DIR, f"qq_PC{i}.png")
    fig_qq, ax_qq = plt.subplots(figsize=(6, 6))
    _draw_qq_ax(ax_qq, arrays, f"QQ plot — PC{i}  (n={n_snps:,})", fontsize=11)
    plt.tight_layout()
    fig_qq.savefig(qq_png, dpi=150)
    plt.close(fig_qq)
    print(f"[PC{i}] QQ plot saved  -> {qq_png}  (λ_GC={lambda_gc:.3f})")

    # -------------------------------------------------------------------
    # Summary row
    # -------------------------------------------------------------------
    summary_rows.append({
        "Trait": f"PC{i}",
        "n_SNPs_tested": n_snps,
        "bonferroni_threshold": bonferroni,
        "n_significant": n_sig,
        "lambda_GC": round(lambda_gc, 4),
        "top_SNP_rs": top_rs,
        "top_SNP_chr": top_chr,
        "top_SNP_pos": top_pos,
        "top_SNP_pval": top_pval,
        "top_SNP_beta": top_beta,
    })

    # -------------------------------------------------------------------
    # Interpretation
    # -------------------------------------------------------------------
    if n_sig > 0:
        sig_chroms = sorted(sig["chr"].astype(int).unique().tolist())
        direction = "positive" if top_beta > 0 else "negative"
        interp_lines.append(
            f"PC{i}: {n_sig} significant SNPs (Bonferroni p<{bonferroni:.2e}) on "
            f"chromosome(s) {sig_chroms}. Top SNP {top_rs} at Chr{top_chr}:{top_pos} "
            f"(p={top_pval:.3e}, beta={top_beta:.4f}, {direction} effect on PC{i} score). "
            f"SNPs clustered in a tight window suggest a QTL; isolated hits may be noise."
        )
    else:
        interp_lines.append(
            f"PC{i}: No SNPs survived Bonferroni correction (threshold p<{bonferroni:.2e}). "
            f"Top nominal hit: {top_rs} at Chr{top_chr}:{top_pos} (p={top_pval:.3e}). "
            f"Weak or no large-effect loci detected for this PC axis."
        )

# ---------------------------------------------------------------------------
# Write summary CSV
# ---------------------------------------------------------------------------
if summary_rows:
    summ_df = pd.DataFrame(summary_rows)
    summ_csv = os.path.join(LMM_DIR, "summary_lmm.csv")
    summ_df.to_csv(summ_csv, index=False)
    print(f"\nSummary CSV -> {summ_csv}")
    print(summ_df.to_string(index=False))

# ---------------------------------------------------------------------------
# Merged QQ panel (2 rows × 4 columns)
# ---------------------------------------------------------------------------
if qq_data:
    fig_panel, axes = plt.subplots(2, 4, figsize=(18, 9))
    axes_flat = axes.flatten()
    for idx, pc in enumerate(sorted(qq_data.keys())):
        n_snps_pc = summary_rows[idx]["n_SNPs_tested"]
        _draw_qq_ax(axes_flat[idx], qq_data[pc],
                    f"PC{pc}  (n={n_snps_pc:,})", fontsize=9)
    fig_panel.suptitle("QQ Plots — LMM GWAS (PC1–PC8)", fontsize=14, y=1.01)
    plt.tight_layout()
    panel_png = os.path.join(LMM_DIR, "qq_panel.png")
    fig_panel.savefig(panel_png, dpi=150, bbox_inches="tight")
    plt.close(fig_panel)
    print(f"\nMerged QQ panel -> {panel_png}")

# ---------------------------------------------------------------------------
# Write interpretation text
# ---------------------------------------------------------------------------
interp_txt = os.path.join(LMM_DIR, "interpretation_lmm.txt")
with open(interp_txt, "w") as fh:
    fh.write("LMM GWAS Interpretation — PC1 to PC8\n")
    fh.write("=" * 60 + "\n\n")
    for line in interp_lines:
        fh.write(line + "\n\n")
print(f"Interpretation -> {interp_txt}")

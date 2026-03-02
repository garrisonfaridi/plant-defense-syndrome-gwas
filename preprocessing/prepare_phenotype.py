"""
prepare_phenotype.py
Align GWAS pseudotraits (PC1–PC8) to the PLINK FAM sample order
and write phenotype.txt for GEMMA.

Run from any directory:
    python preprocessing/prepare_phenotype.py
"""

import os
import pandas as pd

# ---------------------------------------------------------------------------
# Paths  (preprocessing/ → repo root)
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

FAM_PATH    = os.path.join(REPO_ROOT, "genotypes", "genotypes_305.fam")
PSEUDO_PATH = os.path.join(REPO_ROOT, "data", "gwas_pseudotraits_PMM.csv")
OUT_PATH    = os.path.join(REPO_ROOT, "data", "phenotype.txt")

pheno_cols = ["PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8"]

fam = pd.read_csv(
    FAM_PATH,
    sep=r"\s+",
    header=None,
    usecols=[0],
    names=["FamID"],
)

pseudo = pd.read_csv(PSEUDO_PATH, index_col="ecotype_id")

out = pseudo[pheno_cols].reindex(fam["FamID"])
out.to_csv(OUT_PATH, sep="\t", header=False, index=False, na_rep="NA")

n_total   = len(out)
n_non_na  = int(fam["FamID"].isin(pseudo.index).sum())
print(f"phenotype.txt written -> {OUT_PATH}")
print(f"  {n_total} rows x {len(pheno_cols)} columns")
print(f"  Non-NA rows : {n_non_na}")
print(f"  NA rows     : {n_total - n_non_na}")

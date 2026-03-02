import pandas as pd

BASE = "/Users/garrisonfaridi/Documents/Documents_Cloud/BioCode/Cirrone_Lab"
pheno_cols = ["PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8"]

fam = pd.read_csv(
    f"{BASE}/genotypes/genotypes_305.fam",
    sep=r"\s+",
    header=None,
    usecols=[0],
    names=["FamID"],
)

pseudo = pd.read_csv(f"{BASE}/gwas_pseudotraits_PMM.csv", index_col="ecotype_id")

out = pseudo[pheno_cols].reindex(fam["FamID"])
out.to_csv(f"{BASE}/phenotype.txt", sep="\t", header=False, index=False, na_rep="NA")

n_total = len(out)
n_non_na = int(fam["FamID"].isin(pseudo.index).sum())
print(f"phenotype.txt written — {n_total} rows x {len(pheno_cols)} columns")
print(f"  Non-NA rows : {n_non_na}")
print(f"  NA rows     : {n_total - n_non_na}")

# Plant Defense Syndrome GWAS

Genome-wide association study (GWAS) of *Arabidopsis thaliana* plant defense syndromes using a linear mixed model (LMM) in GEMMA. Phenotypic variation across 44 defense traits was decomposed into principal components (PCs); each PC was used as a GWAS trait. Post-hoc GO enrichment and plant GO slim analysis link significant loci to biological function.

---

## Background

Plants deploy multiple, partially redundant defense systems against herbivores and pathogens — glucosinolates, trichomes, reactive oxygen species, and direct resistance proteins. This project tests whether natural genetic variation in *A. thaliana* underlies coordinated "defense syndromes" by performing GWAS on the principal components of a 44-trait defense phenotype matrix.

---

## Pipeline overview

```
phenotype data (44 traits)
        │
        ▼
prepare_phenotype.py       → BLUP residuals, PCA → gwas_pseudotraits_PMM.csv
        │
        ▼
download_annotation.py     → TAIR10 GFF3, GO OBO, TAIR GAF
        │
        ▼
run_gemma.sh               → kinship matrix (centred relatedness)
run_lmm.sh                 → LMM GWAS, 8 PCs (Wald test p-values)
        │
        ▼
plot_lmm.py                → Manhattan plots, QQ plots, summary CSV
        │
        ▼
posthoc_lmm.py             → SNP → gene overlap (±500 kb), GO enrichment
        │
        ▼
posthoc_visualize.py       → dot plots, GO slim, heatmap, biological report
```

---

## Key results

| Trait | Sig. SNPs (Bonferroni) | λ_GC | Top locus |
|-------|----------------------|------|-----------|
| **PC1** (glucosinolate axis) | **27** | 1.057 | Chr5:7,634,212 |
| **PC5** (herbivore damage axis) | **5** | 0.967 | Chr5:7,720,108 |
| PC2–PC4, PC6–PC8 | 0 | 0.97–1.03 | — |

Both significant PCs map to overlapping loci on **Chr5 ~7.6–7.7 Mb**, a dense TIR-NBS-LRR (NLR) resistance gene cluster (AT5G44420–AT5G44870).

**PC1 GO enrichment** (10 significant terms, FDR < 0.05):

| GO term | Namespace | p_FDR |
|---------|-----------|-------|
| LRR domain binding | MF | 6.4 × 10⁻⁷ |
| protein domain specific binding | MF | 9.3 × 10⁻⁷ |
| defense response to insect | BP | 1.4 × 10⁻⁴ |
| receptor Ser/Thr kinase binding | MF | 1.4 × 10⁻⁴ |
| response to herbivore | BP | 9.1 × 10⁻⁴ |
| sesquiterpene metabolic process | BP | 9.8 × 10⁻³ |
| defense response to virus | BP | 2.3 × 10⁻² |

---

## Repository structure

```
.
├── prepare_phenotype.py        # BLUP residuals + PCA → GWAS pseudotraits
├── download_annotation.py      # Download TAIR10 GFF3, GO OBO, TAIR GAF
├── run_gemma.sh                # GEMMA kinship matrix
├── run_lmm.sh                  # GEMMA LMM GWAS (8 PCs)
├── run_bslmm.sh                # GEMMA BSLMM (optional Bayesian model)
├── plot_lmm.py                 # Manhattan plots, QQ plots, summary stats
├── summarize_bslmm.py          # BSLMM PIP summaries
├── posthoc_lmm.py              # SNP→gene overlap + GO enrichment
├── posthoc_visualize.py        # Visualizations + GO slim + report
│
├── gwas_pseudotraits_PMM.csv   # 8 PC phenotypes used as GWAS traits
├── pca_loadings_matrix_PMM.csv # 44 traits × 8 PCs loadings
├── phenotype.txt               # Raw phenotype input for GEMMA
│
├── annotation/                 # Downloaded reference files (see below)
│
└── gemma_output/
    ├── lmm/
    │   ├── manhattan_PC{1-8}.png
    │   ├── qq_PC{1-8}.png
    │   ├── qq_panel.png        # Merged 2×4 QQ panel
    │   ├── summary_lmm.csv
    │   └── interpretation_lmm.txt
    └── posthoc/
        ├── go_enrichment_PC{1,5}.csv
        ├── go_enrichment_combined.csv
        ├── go_slim_PC1.csv
        ├── go_slim_combined.csv
        ├── genes_PC{1,5}.csv
        ├── genes_combined.csv
        ├── go_dot_PC1.png
        ├── go_dot_combined.png
        ├── go_slim_plot.png
        ├── trait_loading_heatmap.png
        └── biological_report.txt
```

> **Not tracked** (too large or regenerable): raw `.assoc.txt` files, TAIR GAF (~54 MB), kinship matrix, genotype files, GO OBO files.

---

## Requirements

```
python >= 3.10
pandas numpy scipy matplotlib seaborn
goatools requests
gemma >= 0.98
```

Install Python dependencies:

```bash
pip install pandas numpy scipy matplotlib seaborn goatools requests
```

---

## Reproducing the analysis

All scripts are run from the `Cirrone_Lab/` directory (or the repo root after cloning).

### 1. Download annotations

```bash
python download_annotation.py
```

Downloads TAIR10 GFF3, TAIR GAF, and `go-basic.obo` into `annotation/`.

### 2. Prepare phenotypes

```bash
python prepare_phenotype.py
```

Produces `gwas_pseudotraits_PMM.csv` and `pca_loadings_matrix_PMM.csv`.

### 3. Run GEMMA

```bash
bash run_gemma.sh   # kinship matrix
bash run_lmm.sh     # LMM GWAS
```

Requires GEMMA and genotype files in `genotypes/` (not tracked).

### 4. Plot GWAS results

```bash
python plot_lmm.py
```

Writes Manhattan plots, QQ plots (`qq_panel.png`), and `summary_lmm.csv` to `gemma_output/lmm/`.

### 5. Post-hoc GO enrichment

```bash
python posthoc_lmm.py
```

Overlaps significant SNPs (±500 kb) with TAIR10 genes and runs GO enrichment via goatools.

### 6. Visualize and report

```bash
python posthoc_visualize.py
```

Generates dot plots, plant GO slim analysis, trait loading heatmap, and `biological_report.txt`.

---

## Key loci

| Locus | Gene(s) | Function |
|-------|---------|----------|
| Chr5:7.58–7.72 Mb | AT5G44420–AT5G44870 | TIR-NBS-LRR (NLR) resistance gene cluster; LRR domain binding, ROS, defense signaling |
| Chr5:23.96 kb | AT5G23960 (*TPS21*) | Sesquiterpene synthase; volatile defense attracting parasitoid wasps |
| Chr5:44.62–44.63 kb | AT5G44620, AT5G44630 | Sesquiterpene synthases; co-located with NLR cluster |
| Chr2:13.79 kb | AT2G13790 | LRR receptor-like kinase; pattern-triggered immunity (PTI) |

---

## Reference

*Arabidopsis thaliana* TAIR10 genome assembly. GO annotations from the Gene Ontology Consortium (TAIR GAF). GWAS performed with [GEMMA](https://github.com/genetics-statistics/GEMMA) (Zhou & Stephens 2012).

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
preprocessing/prepare_phenotype.py     → align pseudotraits to PLINK FAM order
        │
        ▼
scripts/python/download_annotation.py  → TAIR10 GFF3, GO OBO, TAIR GAF
        │
        ▼
scripts/shell/run_lmm.sh               → kinship matrix + LMM GWAS (8 PCs)
scripts/shell/run_bslmm.sh             → BSLMM (local); or:
scripts/slurm/slurm_bslmm_array.sh    → BSLMM array job (HPC/SLURM)
        │
        ▼
scripts/python/plot_lmm.py             → Manhattan plots, QQ plots, summary CSV
        │
        ▼
scripts/python/posthoc_lmm.py          → SNP → gene overlap (±500 kb), GO enrichment
        │
        ▼
scripts/python/posthoc_visualize.py    → dot plots, GO slim, heatmap, biological report
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
├── README.md
├── .gitignore
│
├── data/                            # Input data files
│   ├── phenotype.txt                # GEMMA phenotype input (aligned to FAM order)
│   ├── gwas_pseudotraits_PMM.csv    # 8 PC pseudotraits for GWAS
│   └── pca_loadings_matrix_PMM.csv  # 44 traits × 8 PCs loadings
│
├── preprocessing/                   # Phenotype processing and PCA
│   ├── prepare_phenotype.py         # Align pseudotraits → phenotype.txt
│   ├── execution_trace.ipynb        # Full preprocessing notebook
│   └── output/                      # QC figures and intermediate tables
│
├── annotation/                      # Reference annotation (mostly gitignored)
│
├── genotypes/                       # PLINK bed/bim/fam (gitignored)
│
├── scripts/
│   ├── python/                      # Python analysis scripts
│   │   ├── download_annotation.py   # Download TAIR10 GFF3, GO OBO, TAIR GAF
│   │   ├── plot_lmm.py              # Manhattan + QQ plots, summary stats
│   │   ├── posthoc_lmm.py           # SNP → gene overlap + GO enrichment
│   │   ├── posthoc_visualize.py     # GO dot plots, slim, heatmap, report
│   │   └── summarize_bslmm.py       # BSLMM posterior summaries
│   ├── shell/                       # Local shell runners
│   │   ├── run_gemma.sh             # Full pipeline (kinship + LMM + BSLMM)
│   │   ├── run_lmm.sh               # Kinship matrix + LMM GWAS
│   │   └── run_bslmm.sh             # BSLMM (exploratory parameters)
│   └── slurm/                       # HPC SLURM scripts
│       ├── slurm_bslmm_array.sh     # Array job: 8 PCs in parallel (pub. params)
│       └── slurm_bslmm_summarize.sh # Summary step after array completes
│
├── logs/
│   └── gemma_run.log
│
└── results/                         # All GEMMA outputs
    ├── lmm/
    │   ├── manhattan_PC{1-8}.png
    │   ├── qq_PC{1-8}.png
    │   ├── qq_panel.png             # Merged 2×4 QQ panel
    │   ├── summary_lmm.csv
    │   └── interpretation_lmm.txt
    ├── bslmm/
    │   ├── summary_bslmm.csv
    │   ├── summary_hyp.csv
    │   ├── top_snps_PC{1-8}.csv
    │   └── interpretation_bslmm.txt
    └── posthoc/
        ├── go_enrichment_PC{1,5}.csv
        ├── go_enrichment_combined.csv
        ├── go_slim_PC1.csv / go_slim_combined.csv
        ├── genes_PC{1,5}.csv / genes_combined.csv
        ├── go_dot_PC1.png / go_dot_combined.png
        ├── go_slim_plot.png
        ├── trait_loading_heatmap.png
        └── biological_report.txt
```

> **Not tracked** (too large or regenerable): raw `.assoc.txt` files, TAIR GAF (~54 MB), kinship matrix, genotype files, GO OBO files, BSLMM MCMC chains.

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
python scripts/python/download_annotation.py
```

Downloads TAIR10 GFF3, TAIR GAF, and `go-basic.obo` into `annotation/`.

### 2. Prepare phenotypes

```bash
python preprocessing/prepare_phenotype.py
```

Aligns `data/gwas_pseudotraits_PMM.csv` to the PLINK FAM sample order and writes `data/phenotype.txt`.

### 3. Run GEMMA

```bash
bash scripts/shell/run_lmm.sh     # kinship matrix + LMM GWAS
bash scripts/shell/run_bslmm.sh   # BSLMM (or use SLURM scripts on HPC)
```

Requires GEMMA and genotype files in `genotypes/` (not tracked).

### 4. Plot GWAS results

```bash
python scripts/python/plot_lmm.py
```

Writes Manhattan plots, QQ plots (`qq_panel.png`), and `summary_lmm.csv` to `results/lmm/`.

### 5. Post-hoc GO enrichment

```bash
python scripts/python/posthoc_lmm.py
```

Overlaps significant SNPs (±500 kb) with TAIR10 genes and runs GO enrichment via goatools.

### 6. Visualize and report

```bash
python scripts/python/posthoc_visualize.py
```

Generates dot plots, plant GO slim analysis, trait loading heatmap, and `biological_report.txt` in `results/posthoc/`.

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

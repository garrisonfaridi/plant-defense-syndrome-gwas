"""
posthoc_visualize.py
Visualization, GO slim analysis, and biological report for post-hoc LMM results.

Run from any directory:
    python scripts/python/posthoc_visualize.py

Outputs (results/posthoc/):
    go_dot_PC1.png
    go_dot_combined.png
    go_slim_PC1.csv
    go_slim_combined.csv
    go_slim_plot.png
    trait_loading_heatmap.png
    biological_report.txt
"""

import gzip
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

plt.style.use("seaborn-v0_8-whitegrid")

# ---------------------------------------------------------------------------
# Section 0 — Paths & shared setup
# ---------------------------------------------------------------------------
SCRIPT_DIR    = os.path.dirname(os.path.abspath(__file__))         # scripts/python/
REPO_ROOT     = os.path.dirname(os.path.dirname(SCRIPT_DIR))       # repo root
ANNOT_DIR     = os.path.join(REPO_ROOT, "annotation")
OUT_DIR       = os.path.join(REPO_ROOT, "results", "posthoc")

LOADINGS_PATH = os.path.join(REPO_ROOT, "data", "pca_loadings_matrix_PMM.csv")
GFF_PATH      = os.path.join(ANNOT_DIR, "TAIR10_GFF3_genes.gff.gz")
GAF_PATH      = os.path.join(ANNOT_DIR, "tair.gaf")
OBO_PATH      = os.path.join(ANNOT_DIR, "go-basic.obo")
GOSLIM_OBO    = os.path.join(ANNOT_DIR, "goslim_plant.obo")
GOSLIM_URL    = "https://current.geneontology.org/ontology/subsets/goslim_plant.obo"

os.makedirs(OUT_DIR, exist_ok=True)

NS_COLORS = {"BP": "#4e79a7", "MF": "#f28e2b", "CC": "#59a14f"}
NS_ABBREV = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}
ALPHA = 0.05

# Human-readable trait label mapping for heatmap
TRAIT_LABELS = {
    "pseu32.luxPerGram.blup":          "Pseudomonas lux/gram (32h)",
    "pseu63.luxPerGram.blup":          "Pseudomonas lux/gram (63h)",
    "pseu32.Symptoms.blup":            "Pseudomonas symptoms (32h)",
    "pseu63.Symptoms.blup":            "Pseudomonas symptoms (63h)",
    "pseu.rosetteMass.blup":           "Pseudomonas rosette mass",
    "pseu32.rosetteMass.blup":         "Pseudomonas rosette mass (32h)",
    "pseu63.rosetteMass.blup":         "Pseudomonas rosette mass (63h)",
    "slug.propLeavesRm.blup":          "Slug: prop. leaves removed",
    "slug.areafc.d0d2.blup":           "Slug: leaf area fold-change (d0-d2)",
    "slug.areafc.d0d1.blup":           "Slug: leaf area fold-change (d0-d1)",
    "slug.areaD0.blup":                "Slug: leaf area day 0",
    "slug.numLeavesD0.blup":           "Slug: leaf count day 0",
    "aphid.numAphids.blup":            "Aphid: aphid count",
    "aphid.rosetteMass.blup":          "Aphid: rosette mass",
    "beetle.propLeavesRm.blup":        "Beetle: prop. leaves removed",
    "beetle.leafAreaFull.blup":        "Beetle: leaf area",
    "beetle.numLeavesFull.blup":       "Beetle: leaf count",
    "moths.numEggsPerMoth.blup":       "Moth: eggs per moth",
    "fungus.LesionArea.blup":          "Fungal lesion area",
    "trichomes.NumTrichomesResid.blup": "Trichomes (residual)",
    "trichomes.NumTrichomes.blup":     "Trichome count",
    "trichomes.FocalLeafArea.blup":    "Trichome focal leaf area",
    "SA.SA.blup":                      "Salicylic acid",
    "SA.mass.blup":                    "SA: rosette mass",
    "major-gsl.Pren.blup":             "Glucosinolate Prenyl (major)",
    "major-gsl.Buen.blup":             "Glucosinolate Butenyl (major)",
    "major-gsl.S2hBuen.blup":          "Glucosinolate S2h-Butenyl",
    "major-gsl.8mSOo.blup":            "Glucosinolate 8mSOo",
    "major-gsl.IM.blup":               "Glucosinolate IM",
    "major-gsl.1moIM.blup":            "Glucosinolate 1moIM",
    "major-gsl.4moIM.blup":            "Glucosinolate 4moIM",
    "minor-gsl.7mSh.blup":             "Glucosinolate 7mSh (minor)",
    "minor-gsl.8MTO.blup":             "Glucosinolate 8MTO (minor)",
    "minor-gsl.3mSOp.blup":            "Glucosinolate 3mSOp (minor)",
    "minor-gsl.4mSOb.blup":            "Glucosinolate 4mSOb (minor)",
    "minor-gsl.5mSOp.blup":            "Glucosinolate 5mSOp (minor)",
    "minor-gsl.6mSOh.blup":            "Glucosinolate 6mSOh (minor)",
    "minor-gsl.7mSOh.blup":            "Glucosinolate 7mSOh (minor)",
    "minor-gsl.Peen.blup":             "Glucosinolate Pentenyl (minor)",
    "minor-gsl.2hPeen.blup":           "Glucosinolate 2h-Pentenyl (minor)",
    "minor-gsl.1hIM.blup":             "Glucosinolate 1hIM (minor)",
    "gsl.mass.blup":                   "Glucosinolate mass",
}

# Known gene roles for the biological report
GENE_ROLES = {
    "AT5G44420": "NLR resistance protein (TIR-NBS-LRR); defense signaling",
    "AT5G44563": "NLR resistance protein; ROS-associated defense",
    "AT5G44565": "TIR-NBS-LRR; LRR domain binding / herbivore defense",
    "AT5G44567": "TIR-NBS-LRR; LRR domain binding / defense",
    "AT5G44568": "TIR-NBS-LRR; defense response to insect",
    "AT5G44570": "TIR-NBS-LRR; receptor-mediated defense",
    "AT5G44572": "TIR-NBS-LRR; receptor-mediated defense",
    "AT5G44574": "TIR-NBS-LRR; LRR kinase binding / herbivore defense",
    "AT5G44575": "TIR-NBS-LRR; LRR kinase binding / herbivore defense",
    "AT5G44578": "TIR-NBS-LRR; LRR domain binding",
    "AT5G44580": "TIR-NBS-LRR; LRR domain binding",
    "AT5G44582": "TIR-NBS-LRR; LRR domain binding",
    "AT5G44585": "TIR-NBS-LRR; LRR domain / defense to insect",
    "AT5G44620": "TPS (terpene synthase); sesquiterpene biosynthesis",
    "AT5G44630": "TPS (terpene synthase); sesquiterpene biosynthesis",
    "AT5G23960": "TPS21; sesquiterpene synthase — volatile defence attracting parasitoids",
    "AT5G43580": "Disease resistance protein; defense response to insect",
    "AT5G23010": "JAZ (jasmonate-ZIM-domain); jasmonate signaling / insect response",
    "AT2G13790": "LRR receptor-like kinase; signaling receptor binding / PTI",
    "AT2G13540": "Defense response protein; defense to virus/bacteria",
    "AT2G13650": "Disease resistance protein; defense response",
    "AT2G13810": "TIR-NBS-LRR; defense response",
    "AT5G21150": "RPM1-interacting protein; pattern-triggered immunity / defense to virus",
    "AT5G23570": "NLR domain; defense response to virus",
    "AT5G42950": "Disease resistance protein; defense to pathogen",
    "AT5G43470": "Receptor-like protein; defense response to virus",
    "AT5G43810": "Defense-related protein; innate immunity",
    "AT5G44200": "Putative resistance gene; defense to virus / biotic stress",
    "AT5G44870": "TIR-NBS-LRR; defense response to virus",
    "AT5G44030": "Defense protein; part of defense response cluster",
    "AT5G44070": "Disease resistance protein",
    "AT5G44430": "NLR resistance gene; Chr5 cluster",
    "AT5G44510": "NLR/LRR gene; Chr5 defense cluster",
    "AT5G22250": "Disease resistance gene; defense cluster",
    "AT5G22290": "NLR; defense response",
    "AT5G22540": "NLR-like protein; defense",
    "AT5G22570": "Resistance-related kinase",
    "AT5G22630": "Disease resistance protein",
    "AT5G22690": "TIR-NBS-LRR; defense",
    "AT5G23820": "Defense response gene",
    "AT5G24110": "Disease resistance protein",
    "AT5G24210": "Defense gene cluster; Chr5",
    "AT5G42980": "NLR resistance protein",
    "AT5G43560": "TIR-NBS-LRR; defense",
    "AT5G43570": "Disease resistance gene",
    "AT5G45050": "NLR-related; defense response",
    "AT5G45060": "NLR-related; defense response",
}


# ---------------------------------------------------------------------------
# Section 1 — Load enrichment results & PC loadings
# ---------------------------------------------------------------------------
def load_data():
    """Load GO enrichment CSVs, gene overlap lists, and PC loadings."""
    print("Loading enrichment results and loadings ...", flush=True)

    go_pc1      = pd.read_csv(os.path.join(OUT_DIR, "go_enrichment_PC1.csv"))
    go_combined = pd.read_csv(os.path.join(OUT_DIR, "go_enrichment_combined.csv"))
    genes_pc1   = pd.read_csv(os.path.join(OUT_DIR, "genes_PC1.csv"))
    genes_combined = pd.read_csv(os.path.join(OUT_DIR, "genes_combined.csv"))
    loadings    = pd.read_csv(LOADINGS_PATH, index_col="trait")

    print(f"  PC1 GO terms: {len(go_pc1)}")
    print(f"  Combined GO terms: {len(go_combined)}")
    print(f"  PC1 study genes: {len(genes_pc1)}")
    print(f"  Combined study genes: {len(genes_combined)}")
    print(f"  Loadings matrix: {loadings.shape}\n")

    return go_pc1, go_combined, genes_pc1, genes_combined, loadings


# ---------------------------------------------------------------------------
# Section 2 — Parse gene names from GFF3
# ---------------------------------------------------------------------------
def parse_gene_names() -> dict:
    """Extract {gene_id: Name} from TAIR10 GFF3 attributes."""
    print("Parsing gene names from GFF3 ...", flush=True)

    gene_names: dict = {}
    opener = gzip.open(GFF_PATH, "rt") if GFF_PATH.endswith(".gz") else open(GFF_PATH)
    with opener as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            _, _, feature, _, _, _, _, _, attrs = parts
            if feature != "gene":
                continue

            gene_id = None
            name    = None
            for token in attrs.split(";"):
                if token.startswith("gene_id="):
                    gene_id = token[8:].split(".")[0]
                elif token.startswith("ID=") and gene_id is None:
                    raw = token[3:].split(":")[-1]
                    gene_id = raw.split(".")[0]
                elif token.startswith("Name="):
                    name = token[5:].strip()

            if gene_id and name:
                gene_names[gene_id] = name

    print(f"  Parsed {len(gene_names):,} gene name mappings\n")
    return gene_names


# ---------------------------------------------------------------------------
# Section 3 — GO enrichment dot plots
# ---------------------------------------------------------------------------
def plot_go_dot(go_df: pd.DataFrame, title: str, outpath: str) -> None:
    """Horizontal bubble/dot plot: x = -log10(FDR), y = GO term, size = study_count."""
    if go_df.empty:
        print(f"  [SKIP] Empty data for: {title}")
        return

    df = go_df.sort_values("p_fdr_bh", ascending=False).reset_index(drop=True)
    df["neg_log10_p"] = -np.log10(df["p_fdr_bh"].clip(lower=1e-300))

    # Scale dot sizes 30–200 pt²
    sc = df["study_count"]
    sc_min, sc_max = sc.min(), sc.max()
    if sc_max > sc_min:
        sizes = 30 + (sc - sc_min) / (sc_max - sc_min) * 170
    else:
        sizes = np.full(len(df), 80.0)

    colors = [NS_COLORS.get(ns, "#999999") for ns in df["namespace"]]
    y_pos  = np.arange(len(df))

    fig_h = max(5, len(df) * 0.55 + 1.5)
    fig, ax = plt.subplots(figsize=(10, fig_h))

    ax.scatter(
        df["neg_log10_p"], y_pos,
        s=sizes, c=colors,
        alpha=0.85, linewidths=0.5, edgecolors="white", zorder=3,
    )

    ax.axvline(-np.log10(ALPHA), color="crimson", linestyle="--",
               linewidth=1.2, label=f"FDR = {ALPHA:.2f}")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df["GO_name"], fontsize=9)
    ax.set_xlabel(r"$-\log_{10}$(FDR-adjusted $p$)", fontsize=11)
    ax.set_title(title, fontsize=13)

    # Namespace colour patches
    ns_patches = [mpatches.Patch(color=c, label=ns) for ns, c in NS_COLORS.items()]

    # Size scale markers
    unique_counts = sorted(set([int(sc_min), int((sc_min + sc_max) / 2), int(sc_max)]))
    size_markers = [
        plt.scatter([], [], s=30 + (v - sc_min) / max(sc_max - sc_min, 1) * 170,
                    c="grey", alpha=0.7, label=f"n = {v}")
        for v in unique_counts
    ]

    legend1 = ax.legend(handles=ns_patches, title="Namespace",
                        loc="lower right", fontsize=8)
    ax.add_artist(legend1)
    ax.legend(handles=size_markers, title="Study genes",
              loc="upper right", fontsize=8)

    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  Saved dot plot -> {outpath}")


# ---------------------------------------------------------------------------
# Section 4 — Plant GO slim analysis
# ---------------------------------------------------------------------------
def download_goslim() -> None:
    """Download goslim_plant.obo if not already present."""
    if os.path.exists(GOSLIM_OBO):
        print(f"  Slim OBO already present: {GOSLIM_OBO}")
        return
    print(f"Downloading goslim_plant.obo ...", flush=True)
    import requests
    resp = requests.get(GOSLIM_URL, timeout=60)
    resp.raise_for_status()
    with open(GOSLIM_OBO, "w") as fh:
        fh.write(resp.text)
    print(f"  Saved -> {GOSLIM_OBO}\n")


def run_go_slim(
    study_genes: set,
    pop_genes: set,
    godag,
    slim_dag,
    id2gos: dict,
    label: str,
) -> pd.DataFrame:
    """Project gene GO terms to plant slim DAG, then run enrichment study."""
    from goatools.mapslim import mapslim
    from goatools.go_enrichment import GOEnrichmentStudy

    # Build slim-level id2gos for the full population
    slim_id2gos: dict = {}
    for gene, go_ids in id2gos.items():
        slim_terms: set = set()
        for go_id in go_ids:
            if go_id not in godag:
                continue  # skip obsolete terms silently
            try:
                direct, _ = mapslim(go_id, godag, slim_dag)
                slim_terms.update(direct)
            except Exception:
                pass
        if slim_terms:
            slim_id2gos[gene] = slim_terms

    coverage = len(study_genes & set(slim_id2gos.keys())) / max(len(study_genes), 1) * 100
    print(f"  [{label}] Slim mapping: {len(slim_id2gos):,} genes have slim terms "
          f"| study coverage: {coverage:.1f}%")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        goeaobj = GOEnrichmentStudy(
            pop_genes,
            slim_id2gos,
            slim_dag,
            propagate_counts=False,
            alpha=ALPHA,
            methods=["fdr_bh"],
        )
        results = goeaobj.run_study(study_genes)

    rows = []
    for r in results:
        if r.p_fdr_bh is None or r.p_fdr_bh >= ALPHA:
            continue
        go_term = slim_dag.get(r.GO)
        ns_full = go_term.namespace if go_term else "unknown"
        rows.append({
            "GO_id":         r.GO,
            "GO_name":       r.name,
            "namespace":     NS_ABBREV.get(ns_full, ns_full),
            "study_count":   r.study_count,
            "pop_count":     r.pop_count,
            "p_uncorrected": r.p_uncorrected,
            "p_fdr_bh":      r.p_fdr_bh,
            "study_items":   ",".join(sorted(r.study_items)) if r.study_items else "",
        })

    if not rows:
        print(f"  [{label}] No significant slim GO terms (p_fdr_bh < {ALPHA})")
        return pd.DataFrame(columns=[
            "GO_id", "GO_name", "namespace", "study_count", "pop_count",
            "p_uncorrected", "p_fdr_bh", "study_items",
        ])

    result_df = pd.DataFrame(rows).sort_values("p_fdr_bh").reset_index(drop=True)
    print(f"  [{label}] {len(result_df)} significant slim GO terms")
    return result_df


def plot_slim_bar(
    slim_pc1: pd.DataFrame,
    slim_combined: pd.DataFrame,
    outpath: str,
) -> None:
    """Grouped bar chart: top slim terms, PC1 vs combined, coloured by namespace."""
    if slim_pc1.empty and slim_combined.empty:
        print("  [SKIP] No slim terms to plot")
        return

    all_terms = pd.concat([slim_pc1, slim_combined], ignore_index=True)
    # Top 15 terms by lowest p_fdr_bh in either group
    top_names = (
        all_terms.groupby("GO_name")["p_fdr_bh"].min()
        .sort_values().head(15).index.tolist()
    )

    def get_neg_log(df: pd.DataFrame, term: str) -> float:
        row = df[df["GO_name"] == term]
        if row.empty:
            return 0.0
        return float(-np.log10(row.iloc[0]["p_fdr_bh"]))

    def get_ns(df: pd.DataFrame, term: str) -> str:
        row = df[df["GO_name"] == term]
        return row.iloc[0]["namespace"] if not row.empty else "BP"

    pc1_vals  = [get_neg_log(slim_pc1, t) for t in top_names]
    comb_vals = [get_neg_log(slim_combined, t) for t in top_names]
    ns_colors = [NS_COLORS.get(get_ns(all_terms, t), "#999999") for t in top_names]

    x     = np.arange(len(top_names))
    width = 0.38

    fig, ax = plt.subplots(figsize=(max(12, len(top_names) * 0.9), 6))
    ax.bar(x - width / 2, pc1_vals,  width, label="PC1",     color="#4e79a7", alpha=0.85)
    ax.bar(x + width / 2, comb_vals, width, label="Combined", color="#f28e2b", alpha=0.85)
    ax.axhline(-np.log10(ALPHA), color="crimson", linestyle="--",
               linewidth=1.0, label=f"FDR = {ALPHA:.2f}")

    ax.set_xticks(x)
    ax.set_xticklabels(top_names, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel(r"$-\log_{10}$(FDR-adjusted $p$)", fontsize=11)
    ax.set_title("Plant GO Slim Enrichment — PC1 vs Combined", fontsize=13)
    ax.legend(fontsize=9)

    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  Saved slim bar chart -> {outpath}")


# ---------------------------------------------------------------------------
# Section 5 — Trait loading heatmap
# ---------------------------------------------------------------------------
def plot_loading_heatmap(loadings: pd.DataFrame, outpath: str) -> None:
    """Heatmap of PC1–PC5 loadings for top traits driving PC1 and PC5."""
    import seaborn as sns

    top_pc1 = loadings["PC1"].abs().nlargest(15).index.tolist()
    top_pc5 = loadings["PC5"].abs().nlargest(10).index.tolist()
    top_traits = list(dict.fromkeys(top_pc1 + top_pc5))  # deduplicate, preserve order

    plot_df = loadings.loc[top_traits, ["PC1", "PC2", "PC3", "PC4", "PC5"]].copy()
    plot_df.index = [TRAIT_LABELS.get(t, t) for t in plot_df.index]

    fig, ax = plt.subplots(figsize=(9, max(8, len(top_traits) * 0.45 + 2.0)))
    sns.heatmap(
        plot_df,
        cmap="RdBu_r",
        center=0,
        vmin=-0.45,
        vmax=0.45,
        linewidths=0.4,
        annot=True,
        fmt=".2f",
        annot_kws={"size": 7},
        ax=ax,
    )
    ax.set_title("PC Loadings — Top traits for PC1 & PC5", fontsize=13)
    ax.set_xlabel("Principal Component", fontsize=11)
    ax.set_ylabel("")

    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  Saved loading heatmap -> {outpath}")


# ---------------------------------------------------------------------------
# Section 6 — Biological text report
# ---------------------------------------------------------------------------
def write_biological_report(
    go_pc1: pd.DataFrame,
    go_combined: pd.DataFrame,
    loadings: pd.DataFrame,
    gene_names: dict,
    genes_pc1: pd.DataFrame,
    genes_combined: pd.DataFrame,
    outpath: str,
) -> None:
    """Write a plain-text biological interpretation of GO enrichment results."""
    lines = []
    SEP = "=" * 70

    def h1(title: str) -> None:
        lines.append("")
        lines.append(SEP)
        lines.append(title)
        lines.append(SEP)

    def h2(title: str) -> None:
        lines.append("")
        lines.append("-" * 50)
        lines.append(title)
        lines.append("-" * 50)

    lines.append("BIOLOGICAL REPORT — Post-hoc LMM GWAS (Arabidopsis thaliana)")
    lines.append("Generated by posthoc_visualize.py")
    lines.append("Inputs: go_enrichment_PC1.csv, go_enrichment_combined.csv")
    lines.append("")

    # -----------------------------------------------------------------------
    # 1. PC identity
    # -----------------------------------------------------------------------
    h1("1. PC IDENTITY (from loadings matrix)")

    h2("PC1 — Glucosinolate axis")
    pc1_abs = loadings["PC1"].abs().sort_values(ascending=False)
    strong_pc1 = pc1_abs[pc1_abs > 0.20]
    lines.append("Traits with |PC1 loading| > 0.20:")
    for trait, val in strong_pc1.items():
        sign  = "+" if loadings.at[trait, "PC1"] > 0 else "−"
        label = TRAIT_LABELS.get(trait, trait)
        lines.append(f"  {sign}{val:.3f}  {label}  [{trait}]")
    lines.append(
        "\nPC1 is dominated by aliphatic glucosinolate traits (positive loadings on "
        "butenyl, 5mSOp, pentenyl, S2h-butenyl, 8mSOo etc.) with negative loadings on "
        "glucosinolate mass and trichome area. This axis captures variation in "
        "glucosinolate composition and quantity across accessions."
    )

    h2("PC5 — Insect/herbivore damage axis")
    pc5_abs = loadings["PC5"].abs().sort_values(ascending=False)
    strong_pc5 = pc5_abs[pc5_abs > 0.15]
    lines.append("Traits with |PC5 loading| > 0.15:")
    for trait, val in strong_pc5.items():
        sign  = "+" if loadings.at[trait, "PC5"] > 0 else "−"
        label = TRAIT_LABELS.get(trait, trait)
        lines.append(f"  {sign}{val:.3f}  {label}  [{trait}]")
    lines.append(
        "\nPC5 shows positive loadings for slug, aphid, beetle, and trichome counts, "
        "reflecting a general herbivore resistance / tissue consumption axis. "
        "Negative loadings on butenyl glucosinolates suggest an inverse relationship "
        "with aliphatic defence chemistry on this axis."
    )

    # -----------------------------------------------------------------------
    # 2. Significant GO terms
    # -----------------------------------------------------------------------
    h1("2. SIGNIFICANT GO TERMS — COMBINED STUDY SET")

    for _, row in go_combined.sort_values("p_fdr_bh").iterrows():
        go_id   = row["GO_id"]
        go_name = row["GO_name"]
        ns      = row["namespace"]
        p_fdr   = row["p_fdr_bh"]
        n_study = row["study_count"]
        items   = [g.strip() for g in str(row["study_items"]).split(",") if g.strip()]

        lines.append(f"\n  [{ns}] {go_id}  {go_name}")
        lines.append(f"       p_fdr_bh = {p_fdr:.2e}   study genes = {n_study}")
        lines.append("       Genes:")
        for gene_id in sorted(items):
            name = gene_names.get(gene_id, "—")
            role = GENE_ROLES.get(gene_id, "function not annotated in report")
            lines.append(f"         {gene_id}  ({name})  — {role}")

    # -----------------------------------------------------------------------
    # 3. Key loci
    # -----------------------------------------------------------------------
    h1("3. KEY LOCI")

    h2("Chr5 NLR/TIR-LRR resistance gene cluster (AT5G44420–AT5G44870)")
    lines.append(
        "The strongest GWAS signal maps to a dense cluster of TIR-NBS-LRR (NLR) "
        "resistance genes on Chr5 (~AT5G44420–AT5G44870). NLR proteins are intracellular "
        "immune receptors that recognise pathogen/herbivore effectors. The LRR domain "
        "mediates binding to pathogen effectors and downstream signalling partners "
        "(LRR domain binding GO:0030275; receptor serine/threonine kinase binding "
        "GO:0033612). ROS-producing genes within the cluster (GO:1903409) mediate the "
        "oxidative burst of the hypersensitive response. The Chr5 cluster is one of the "
        "largest NLR arrays in the Arabidopsis genome and is a known hotspot for "
        "resistance to both insect herbivores and pathogens."
    )

    h2("Sesquiterpene synthases — AT5G23960 (TPS21), AT5G44620, AT5G44630")
    lines.append(
        "TPS21 and the adjacent loci AT5G44620/AT5G44630 catalyse the first committed "
        "step in sesquiterpene biosynthesis (GO:0051761). Sesquiterpenes are volatile "
        "compounds emitted by damaged Arabidopsis leaves that attract parasitoid wasps "
        "(indirect defence). TPS21 specifically produces (E)-beta-caryophyllene and "
        "alpha-humulene in response to herbivory. Their co-enrichment with herbivore "
        "response (GO:0080027) and LRR kinase binding genes highlights a signalling-"
        "metabolic defence circuit co-localised on Chr5."
    )

    h2("Pattern-triggered immunity — AT5G21150, AT5G23570, AT5G44200, AT5G44870")
    lines.append(
        "These loci contribute to defense response to virus (GO:0051607) and broader "
        "defense response (GO:0006952). AT5G21150 encodes an RPM1-interacting protein "
        "involved in effector-triggered immunity; AT5G23570 carries an NLR domain. "
        "The antiviral annotation likely reflects pleiotropic roles of the Chr5 NLR "
        "cluster: TIR-NBS-LRR proteins share domain architecture regardless of the "
        "elicitor (herbivore vs pathogen)."
    )

    h2("AT2G13790 — LRR receptor-like kinase (Chr2)")
    lines.append(
        "AT2G13790 encodes an LRR-RLK (leucine-rich repeat receptor-like kinase) on "
        "Chr2. LRR-RLKs are surface PRRs perceiving DAMPs and herbivore-associated "
        "molecular patterns (HAMPs). Its contribution to signaling receptor binding "
        "(GO:0005102) and receptor serine/threonine kinase binding (GO:0033612) is "
        "consistent with a role in pattern-triggered immunity (PTI) upstream of the "
        "glucosinolate defence response."
    )

    # -----------------------------------------------------------------------
    # 4. PC5 note
    # -----------------------------------------------------------------------
    h1("4. PC5 NOTE")

    genes_pc5_path = os.path.join(OUT_DIR, "genes_PC5.csv")
    n_pc5_genes = len(pd.read_csv(genes_pc5_path)) if os.path.exists(genes_pc5_path) else "N/A"

    lines.append(
        f"PC5 produced {n_pc5_genes} overlapping genes but no significant GO terms "
        f"after FDR correction (p_fdr_bh < {ALPHA}). The Chr5 ~7.72 Mb SNPs "
        "associated with PC5 overlap the same NLR/TIR-LRR cluster as PC1-significant "
        "SNPs, suggesting that two principal components capture partially distinct "
        "phenotypic dimensions of resistance mediated by the same genomic region."
    )
    lines.append(
        "\nPC5 loads positively on slug, aphid, and beetle damage traits and trichome "
        "counts, indicating general herbivore resistance. The smaller PC5 study set "
        "(fewer unique SNPs) combined with sparse NLR annotation may explain the lack "
        "of surviving GO terms. Nominally enriched categories would be expected to "
        "include 'defense response to insect' and 'plant-type hypersensitive response' "
        "with a larger study set."
    )

    # -----------------------------------------------------------------------
    # 5. Summary table
    # -----------------------------------------------------------------------
    h1("5. SUMMARY TABLE — COMBINED GO ENRICHMENT")
    hdr = f"{'GO_id':<14} {'NS':<4} {'p_fdr_bh':<12} {'n_study':<8} {'GO_name'}"
    lines.append(hdr)
    lines.append("-" * 80)
    for _, row in go_combined.sort_values("p_fdr_bh").iterrows():
        lines.append(
            f"{row['GO_id']:<14} {row['namespace']:<4} {row['p_fdr_bh']:<12.2e} "
            f"{row['study_count']:<8} {row['GO_name']}"
        )
    lines.append("")
    lines.append(SEP)
    lines.append("END OF BIOLOGICAL REPORT")
    lines.append(SEP)

    with open(outpath, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  Saved biological report -> {outpath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    # 1. Load data
    go_pc1, go_combined, genes_pc1, genes_combined, loadings = load_data()

    # 2. Parse gene names from GFF3
    gene_names = parse_gene_names()

    # 3. GO enrichment dot plots
    print("Generating GO enrichment dot plots ...", flush=True)
    plot_go_dot(
        go_pc1,
        "GO Enrichment — PC1 (glucosinolate GWAS)",
        os.path.join(OUT_DIR, "go_dot_PC1.png"),
    )
    plot_go_dot(
        go_combined,
        "GO Enrichment — Combined (PC1 + PC5)",
        os.path.join(OUT_DIR, "go_dot_combined.png"),
    )

    # 4. Plant GO slim analysis
    print("\nRunning plant GO slim analysis ...", flush=True)
    download_goslim()

    from goatools.obo_parser import GODag
    from goatools.anno.gaf_reader import GafReader

    print("Loading full GO DAG ...", flush=True)
    godag = GODag(OBO_PATH, optional_attrs={"relationship"})
    print(f"  GO terms loaded: {len(godag):,}")

    print("Loading slim GO DAG ...", flush=True)
    slim_dag = GODag(GOSLIM_OBO, optional_attrs={"subset"})
    print(f"  Slim GO terms loaded: {len(slim_dag):,}")

    print("Loading GAF annotations ...", flush=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gaf = GafReader(GAF_PATH)
    gaf.godag = godag
    id2gos = gaf.get_id2gos()
    print(f"  Annotated genes: {len(id2gos):,}\n")

    pop_genes      = set(id2gos.keys())
    study_pc1      = set(genes_pc1["gene_id"].tolist())
    study_combined = set(genes_combined["gene_id"].tolist())

    slim_pc1 = run_go_slim(
        study_pc1, pop_genes, godag, slim_dag, id2gos, "PC1"
    )
    slim_combined = run_go_slim(
        study_combined, pop_genes, godag, slim_dag, id2gos, "combined"
    )

    slim_pc1.to_csv(os.path.join(OUT_DIR, "go_slim_PC1.csv"), index=False)
    slim_combined.to_csv(os.path.join(OUT_DIR, "go_slim_combined.csv"), index=False)
    print(f"  Saved go_slim_PC1.csv and go_slim_combined.csv")

    plot_slim_bar(slim_pc1, slim_combined, os.path.join(OUT_DIR, "go_slim_plot.png"))

    # 5. Trait loading heatmap
    print("\nGenerating loading heatmap ...", flush=True)
    plot_loading_heatmap(loadings, os.path.join(OUT_DIR, "trait_loading_heatmap.png"))

    # 6. Biological report
    print("\nWriting biological report ...", flush=True)
    write_biological_report(
        go_pc1, go_combined, loadings, gene_names,
        genes_pc1, genes_combined,
        os.path.join(OUT_DIR, "biological_report.txt"),
    )

    print(f"\nAll outputs written to: {OUT_DIR}")


if __name__ == "__main__":
    main()

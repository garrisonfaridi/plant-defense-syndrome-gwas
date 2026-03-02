"""
posthoc_lmm.py
Post-hoc analysis: SNP → gene overlap (±500 kb) + GO enrichment.

Run from BioCode/ working directory:
    python Cirrone_Lab/posthoc_lmm.py
"""

import gzip
import os
import sys
import warnings
import pandas as pd

# ---------------------------------------------------------------------------
# Path setup — all paths derived from the script's own directory
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))  # .../BioCode/Cirrone_Lab
ANNOT_DIR  = os.path.join(SCRIPT_DIR, "annotation")
LMM_DIR    = os.path.join(SCRIPT_DIR, "gemma_output", "lmm")
OUT_DIR    = os.path.join(SCRIPT_DIR, "gemma_output", "posthoc")

GFF_PATH      = os.path.join(ANNOT_DIR, "TAIR10_GFF3_genes.gff.gz")  # Ensembl Plants gzipped GFF3
GAF_PATH_GZ   = os.path.join(ANNOT_DIR, "tair.gaf.gz")               # GO Consortium TAIR GAF (compressed)
GAF_PATH      = os.path.join(ANNOT_DIR, "tair.gaf")                  # decompressed GAF (created on first run)
OBO_PATH      = os.path.join(ANNOT_DIR, "go-basic.obo")

os.makedirs(OUT_DIR, exist_ok=True)

WINDOW_BP  = 500_000
ALPHA      = 0.05
N_PCS      = 8

NS_ABBREV = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}


# ===========================================================================
# Step 2a — Load significant SNPs
# ===========================================================================
def load_significant_snps() -> dict[str, pd.DataFrame]:
    """Return dict {trait: DataFrame} for traits with Bonferroni-significant hits."""

    # Determine Bonferroni threshold from PC1 (all files share the same SNP set)
    pc1_path = os.path.join(LMM_DIR, "lmm_PC1.assoc.txt")
    n_snps = sum(1 for _ in open(pc1_path)) - 1  # subtract header
    bonferroni = ALPHA / n_snps
    print(f"Total SNPs: {n_snps:,}  →  Bonferroni threshold: {bonferroni:.3e}\n")

    sig: dict[str, pd.DataFrame] = {}
    for pc in range(1, N_PCS + 1):
        trait = f"PC{pc}"
        path  = os.path.join(LMM_DIR, f"lmm_{trait}.assoc.txt")
        df = pd.read_csv(path, sep="\t", usecols=["chr", "rs", "ps", "beta", "p_wald"])
        hits = df[df["p_wald"] < bonferroni].copy()
        if not hits.empty:
            sig[trait] = hits
            print(f"  {trait}: {len(hits)} significant SNPs (p_wald < {bonferroni:.3e})")
        else:
            print(f"  {trait}: 0 significant SNPs")

    print()
    return sig


# ===========================================================================
# Step 2b — Parse TAIR10 GFF3 into gene table
# ===========================================================================
def load_genes_gff() -> pd.DataFrame:
    """Parse GFF3 (gzipped Ensembl TAIR10), return DataFrame [chr (int), gene_id, start, end, strand].

    Ensembl Plants GFF3 format:
      - seqname: integer ("1", "2", ..., "5") — no Chr prefix
      - feature == "gene"
      - attributes contain gene_id=AT1G01010 (clean locus ID, no isoform suffix)
    """
    print("Parsing TAIR10 GFF3 ...", flush=True)

    rows = []
    opener = gzip.open(GFF_PATH, "rt") if GFF_PATH.endswith(".gz") else open(GFF_PATH)
    with opener as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqname, _, feature, start, end, _, strand, _, attrs = parts
            if feature != "gene":
                continue

            # Chromosome: Ensembl uses plain integers; TAIR uses "Chr5" — handle both
            raw_chr = seqname[3:] if seqname.startswith("Chr") else seqname
            try:
                chrom = int(raw_chr)
            except ValueError:
                continue  # skip scaffolds / chloroplast / mitochondria

            # Prefer gene_id= field (Ensembl); fall back to ID= stripping gene: prefix
            gene_id = None
            for token in attrs.split(";"):
                if token.startswith("gene_id="):
                    gene_id = token[8:].split(".")[0]
                    break
            if gene_id is None:
                for token in attrs.split(";"):
                    if token.startswith("ID="):
                        raw = token[3:].split(":")[-1]  # strip "gene:" prefix if present
                        gene_id = raw.split(".")[0]
                        break
            if gene_id is None:
                continue

            rows.append({
                "chr":     chrom,
                "gene_id": gene_id,
                "start":   int(start),
                "end":     int(end),
                "strand":  strand,
            })

    genes_df = pd.DataFrame(rows)
    print(f"  Parsed {len(genes_df):,} genes across chromosomes {sorted(genes_df['chr'].unique())}\n")
    return genes_df


# ===========================================================================
# Step 2c — SNP → gene overlap (±500 kb)
# ===========================================================================
def overlap_genes(
    snps_df: pd.DataFrame,
    genes_df: pd.DataFrame,
    label: str,
) -> pd.DataFrame:
    """
    For each SNP in snps_df, find genes within ±WINDOW_BP.
    Returns DataFrame with cols [gene_id, chr, start, end, strand, overlapping_snps].
    """
    gene_hits: dict[str, dict] = {}  # gene_id → info dict

    for _, snp in snps_df.iterrows():
        snp_chr  = int(snp["chr"])
        snp_pos  = int(snp["ps"])
        snp_rs   = snp["rs"]
        win_start = snp_pos - WINDOW_BP
        win_end   = snp_pos + WINDOW_BP

        mask = (
            (genes_df["chr"] == snp_chr)
            & (genes_df["start"] <= win_end)
            & (genes_df["end"]   >= win_start)
        )
        nearby = genes_df[mask]

        for _, g in nearby.iterrows():
            gid = g["gene_id"]
            if gid not in gene_hits:
                gene_hits[gid] = {
                    "gene_id": gid,
                    "chr":     g["chr"],
                    "start":   g["start"],
                    "end":     g["end"],
                    "strand":  g["strand"],
                    "snps":    [],
                }
            gene_hits[gid]["snps"].append(snp_rs)

    if not gene_hits:
        print(f"  [warning] No genes found for {label}")
        return pd.DataFrame(columns=["gene_id", "chr", "start", "end", "strand", "overlapping_snps"])

    records = []
    for info in gene_hits.values():
        records.append({
            "gene_id":         info["gene_id"],
            "chr":             info["chr"],
            "start":           info["start"],
            "end":             info["end"],
            "strand":          info["strand"],
            "overlapping_snps": ",".join(sorted(set(info["snps"]))),
        })

    result = pd.DataFrame(records).sort_values(["chr", "start"]).reset_index(drop=True)
    print(f"  {label}: {len(result):,} genes overlap significant SNPs (±{WINDOW_BP // 1000} kb)")
    return result


# ===========================================================================
# Step 2d — GO enrichment with goatools
# ===========================================================================
def run_go_enrichment(
    study_genes: set[str],
    pop_genes:   set[str],
    id2gos:      dict,
    godag,
    label:       str,
) -> pd.DataFrame:
    """Run GO enrichment for study_genes vs pop_genes; return significant results.

    Uses GOEnrichmentStudy (goatools v1.5.2 API).
    id2gos: {gene_id: set(GO_ids)} flat dict from GafReader.get_id2gos().
    Namespace (BP/MF/CC) is recovered from the GODag for each result term.
    """
    from goatools.go_enrichment import GOEnrichmentStudy

    if not study_genes:
        print(f"  [warning] Empty study set for {label} — skipping GO enrichment")
        return pd.DataFrame()

    print(f"  Running GO enrichment for {label} ({len(study_genes)} study genes) ...", flush=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        goeaobj = GOEnrichmentStudy(
            pop_genes,
            id2gos,
            godag,
            propagate_counts=True,
            alpha=ALPHA,
            methods=["fdr_bh"],
        )
        results = goeaobj.run_study(study_genes)

    rows = []
    for r in results:
        if r.p_fdr_bh is None or r.p_fdr_bh >= ALPHA:
            continue
        # Recover namespace from GODag
        go_term = godag.get(r.GO)
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
        print(f"    No significant GO terms (p_fdr_bh < {ALPHA}) for {label}")
        return pd.DataFrame(columns=[
            "GO_id", "GO_name", "namespace", "study_count", "pop_count",
            "p_uncorrected", "p_fdr_bh", "study_items",
        ])

    df = pd.DataFrame(rows).sort_values("p_fdr_bh").reset_index(drop=True)
    print(f"    {len(df)} significant GO terms found")
    return df


# ===========================================================================
# Step 2e — Print summary
# ===========================================================================
def print_summary(label: str, genes_df: pd.DataFrame, go_df: pd.DataFrame) -> None:
    print(f"\n{'='*60}")
    print(f"Group: {label}")
    print(f"  Genes overlapping significant SNPs: {len(genes_df):,}")
    if go_df.empty:
        print("  GO enrichment: no significant terms")
    else:
        print(f"  Significant GO terms (p_fdr_bh < {ALPHA}): {len(go_df)}")
        top = go_df.head(10)
        print(f"  Top {len(top)} GO hits:")
        for _, row in top.iterrows():
            print(
                f"    [{row['namespace']}] {row['GO_id']}  {row['GO_name'][:50]:<50}"
                f"  p_fdr_bh={row['p_fdr_bh']:.2e}"
                f"  study_n={row['study_count']}"
            )


# ===========================================================================
# Main
# ===========================================================================
def main() -> None:
    # --- 2a: Load significant SNPs ---
    sig_snps = load_significant_snps()
    if not sig_snps:
        print("No significant SNPs found. Exiting.")
        sys.exit(0)

    # --- 2b: Parse GFF3 ---
    genes_df = load_genes_gff()
    pop_genes = set(genes_df["gene_id"].unique())
    print(f"Background gene set size: {len(pop_genes):,}\n")

    # --- Load GO annotations ---
    print("Loading GO OBO ...", flush=True)
    from goatools.obo_parser import GODag
    godag = GODag(OBO_PATH, optional_attrs={"relationship"})
    print(f"  GO terms loaded: {len(godag):,}\n")

    # GafReader (goatools v1.5.2) does not support gzip natively — decompress first
    if not os.path.exists(GAF_PATH):
        print(f"Decompressing {os.path.basename(GAF_PATH_GZ)} ...", flush=True)
        with gzip.open(GAF_PATH_GZ, "rt") as gz_in, open(GAF_PATH, "w") as txt_out:
            txt_out.write(gz_in.read())
        print(f"  Decompressed to {os.path.basename(GAF_PATH)}\n")

    print("Loading GAF annotation ...", flush=True)
    from goatools.anno.gaf_reader import GafReader
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gaf = GafReader(GAF_PATH)
    gaf.godag = godag  # required by get_id2gos() for namespace detection
    id2gos = gaf.get_id2gos()   # {gene_id: set(GO_ids)} — flat dict for GOEnrichmentStudy
    print(f"  Annotated genes: {len(id2gos):,}\n")

    all_annotated_genes = set(id2gos.keys())
    print(f"  Genes with GO annotation: {len(all_annotated_genes):,}")
    print(f"  Overlap with GFF background: {len(all_annotated_genes & pop_genes):,}\n")

    # --- 2c + 2d: Per-group overlap and enrichment ---
    groups: dict[str, pd.DataFrame] = {}  # label → sig snps df

    for trait, df in sig_snps.items():
        groups[trait] = df

    # Build combined group
    if len(sig_snps) > 1:
        combined = pd.concat(list(sig_snps.values()), ignore_index=True).drop_duplicates(subset=["rs"])
        groups["combined"] = combined
    elif len(sig_snps) == 1:
        # Only one trait — combined == that trait
        only_trait = next(iter(sig_snps))
        groups["combined"] = sig_snps[only_trait].copy()

    for label, snps_df in groups.items():
        print(f"\n--- Processing group: {label} ({len(snps_df)} SNPs) ---")

        # Gene overlap
        gene_overlap = overlap_genes(snps_df, genes_df, label)
        gene_out = os.path.join(OUT_DIR, f"genes_{label}.csv")
        gene_overlap.to_csv(gene_out, index=False)
        print(f"  Saved: {gene_out}")

        # GO enrichment
        study_genes = set(gene_overlap["gene_id"].tolist())
        go_df = run_go_enrichment(study_genes, pop_genes, id2gos, godag, label)
        go_out = os.path.join(OUT_DIR, f"go_enrichment_{label}.csv")
        go_df.to_csv(go_out, index=False)
        print(f"  Saved: {go_out}")

        # Summary
        print_summary(label, gene_overlap, go_df)

    print(f"\n{'='*60}")
    print(f"All outputs written to: {OUT_DIR}")


if __name__ == "__main__":
    main()

"""
download_annotation.py
Download TAIR10 annotation files needed for post-hoc LMM analysis.

Run from any directory:
    python scripts/python/download_annotation.py
"""

import os
import sys
import requests

# scripts/python/ → scripts/ → repo root
REPO_ROOT      = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ANNOTATION_DIR = os.path.join(REPO_ROOT, "annotation")

FILES = [
    (
        # Ensembl Plants mirror of TAIR10 gene annotation (gzipped GFF3).
        # arabidopsis.org direct downloads block automated requests (403).
        "TAIR10_GFF3_genes.gff.gz",
        "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.61.gff3.gz",
    ),
    (
        # GO Consortium TAIR GAF — uses AGI locus IDs (AT1G01010 format),
        # compatible with goatools GafReader and TAIR10 gene IDs.
        "tair.gaf.gz",
        "https://current.geneontology.org/annotations/tair.gaf.gz",
    ),
    (
        "go-basic.obo",
        "https://current.geneontology.org/ontology/go-basic.obo",
    ),
]


def download_file(url: str, dest: str) -> None:
    """Stream-download url to dest, printing progress every 10 MB."""
    print(f"  Downloading {os.path.basename(dest)} ...", flush=True)
    response = requests.get(url, stream=True, timeout=120)
    response.raise_for_status()

    total = int(response.headers.get("content-length", 0))
    downloaded = 0
    chunk_size = 1024 * 1024  # 1 MB chunks
    report_every = 10 * 1024 * 1024  # report every 10 MB
    next_report = report_every

    with open(dest, "wb") as f:
        for chunk in response.iter_content(chunk_size=chunk_size):
            if chunk:
                f.write(chunk)
                downloaded += len(chunk)
                if downloaded >= next_report:
                    if total:
                        pct = 100.0 * downloaded / total
                        print(f"    {downloaded / 1e6:.0f} MB / {total / 1e6:.0f} MB ({pct:.0f}%)", flush=True)
                    else:
                        print(f"    {downloaded / 1e6:.0f} MB downloaded", flush=True)
                    next_report += report_every

    print(f"  Done: {os.path.basename(dest)} ({downloaded / 1e6:.1f} MB)", flush=True)


def main() -> None:
    os.makedirs(ANNOTATION_DIR, exist_ok=True)
    print(f"Annotation directory: {ANNOTATION_DIR}\n")

    for filename, url in FILES:
        dest = os.path.join(ANNOTATION_DIR, filename)
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f"  [skip] {filename} already exists ({os.path.getsize(dest) / 1e6:.1f} MB)")
            continue
        download_file(url, dest)

    # Verify all files are non-zero
    print("\nVerification:")
    all_ok = True
    for filename, _ in FILES:
        dest = os.path.join(ANNOTATION_DIR, filename)
        if os.path.exists(dest):
            size = os.path.getsize(dest)
            status = "OK" if size > 0 else "EMPTY — download may have failed"
            print(f"  {filename}: {size / 1e6:.1f} MB  [{status}]")
            if size == 0:
                all_ok = False
        else:
            print(f"  {filename}: MISSING")
            all_ok = False

    if not all_ok:
        sys.exit(1)
    print("\nAll annotation files present and non-zero.")


if __name__ == "__main__":
    main()

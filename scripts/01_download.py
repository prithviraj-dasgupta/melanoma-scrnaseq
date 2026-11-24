"""Download GSE72056 from GEO and the GENCODE v19 gene annotation, recording size and checksum."""

from pathlib import Path

from melscrna.download import GENCODE_GTF_URL, GEO_URL, download, md5sum

ROOT = Path(__file__).resolve().parents[1]
FILES = {
    GEO_URL: ROOT / "data" / "raw" / "GSE72056_melanoma_single_cell_revised_v2.txt.gz",
    GENCODE_GTF_URL: ROOT / "data" / "reference" / "gencode.v19.annotation.gtf.gz",
}

if __name__ == "__main__":
    for url, dest in FILES.items():
        path = download(url, dest)
        print(f"File: {path.relative_to(ROOT)}")
        print(f"Size: {path.stat().st_size / 1e6:.1f} MB")
        print(f"MD5:  {md5sum(path)}\n")

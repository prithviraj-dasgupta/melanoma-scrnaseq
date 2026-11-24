"""Phase 1: download GSE72056 from GEO and record its size and checksum."""

from pathlib import Path

from melscrna.download import GEO_URL, download, md5sum

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "data" / "raw" / "GSE72056_melanoma_single_cell_revised_v2.txt.gz"

if __name__ == "__main__":
    path = download(GEO_URL, DEST)
    size_mb = path.stat().st_size / 1e6
    print(f"File: {path.relative_to(ROOT)}")
    print(f"Size: {size_mb:.1f} MB")
    print(f"MD5:  {md5sum(path)}")

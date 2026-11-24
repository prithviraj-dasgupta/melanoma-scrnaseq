"""Scripted download of the raw data, so nothing has to be fetched by hand."""

from __future__ import annotations

import hashlib
from pathlib import Path

import requests

GEO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72056/suppl/"
    "GSE72056_melanoma_single_cell_revised_v2.txt.gz"
)

# Gene coordinates for copy-number inference. GENCODE v19 is on hg19/GRCh37, the genome build
# of the 2016 data, so gene symbols match as closely as possible.
GENCODE_GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/"
    "gencode.v19.annotation.gtf.gz"
)


def md5sum(path: Path, chunk_size: int = 1 << 20) -> str:
    """Return the MD5 checksum of a file, read in chunks so memory use stays small."""
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, dest: Path, overwrite: bool = False, timeout: int = 60) -> Path:
    """Download ``url`` to ``dest``.

    Skips the download if ``dest`` already exists (unless ``overwrite``). Writes to a temporary
    ``.part`` file first, so an interrupted download never leaves a truncated file that looks
    complete.
    """
    dest = Path(dest)
    if dest.exists() and not overwrite:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    with requests.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with open(tmp, "wb") as handle:
            for chunk in response.iter_content(chunk_size=1 << 20):
                handle.write(chunk)
    tmp.replace(dest)
    return dest

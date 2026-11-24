import gzip
from pathlib import Path

import pytest

HEADER = "Cell\tc1\tc2\tc3\tc4\n"
META = (
    "tumor\t79\t79\t80\t80\n"
    "malignant(1=no,2=yes,0=unresolved)\t2\t1\t1\t0\n"
    "non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAFs;6=NK)\t0\t1\t3\t0\n"
)
GENES = (
    "MITF\t8.1\t0\t0\t3.2\n"
    "CD3E\t0\t7.5\t0\t0\n"
    "CD68\t0\t0\t6.9\t0\n"
    "MT-CO1\t9.0\t8.0\t8.5\t7.0\n"
    "MITF\t1.0\t0\t0\t0\n"
)


@pytest.fixture
def tiny_tirosh(tmp_path: Path) -> Path:
    """A 5-gene x 4-cell file in the same layout as GSE72056 (one duplicated gene)."""
    path = tmp_path / "tiny.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(HEADER + META + GENES)
    return path

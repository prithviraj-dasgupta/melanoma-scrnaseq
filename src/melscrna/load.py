"""Read the Tirosh et al. 2016 melanoma matrix (GEO GSE72056) into an AnnData object.

File layout (one tab-separated table, genes x cells):

- row 1: column headers, i.e. cell IDs
- row 2: ``tumor``: patient / tumour number
- row 3: ``malignant(1=no,2=yes,0=unresolved)``
- row 4: ``non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAFs;6=NK)``; 0 for
  malignant and unresolved cells
- remaining rows: one gene per row, values log2(TPM/10 + 1)
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

N_META_ROWS = 3

MALIGNANT_CODES = {0: "unresolved", 1: "non-malignant", 2: "malignant"}
CELL_TYPE_CODES = {
    0: "none",
    1: "T cell",
    2: "B cell",
    3: "Macrophage",
    4: "Endothelial",
    5: "CAF",
    6: "NK cell",
}


def _decode(values: pd.Series, codes: dict[int, str], name: str) -> pd.Categorical:
    as_int = values.astype(float).round().astype(int)
    unknown = sorted(set(as_int) - set(codes))
    if unknown:
        raise ValueError(f"Unexpected codes in '{name}' row: {unknown}")
    return pd.Categorical(as_int.map(codes), categories=list(codes.values()))


def read_tirosh(path: str | Path) -> ad.AnnData:
    """Return cells x genes AnnData with author labels in ``.obs``.

    ``.X`` holds the values exactly as deposited (log2(TPM/10 + 1)); nothing is normalised
    here. Duplicate gene symbols are made unique (``GENE``, ``GENE-1``) and counted in
    ``.uns["n_duplicate_genes"]``.
    """
    table = pd.read_csv(path, sep="\t", index_col=0)
    if table.shape[0] <= N_META_ROWS:
        raise ValueError("File has no gene rows after the 3 metadata rows")

    meta = table.iloc[:N_META_ROWS]
    expr = table.iloc[N_META_ROWS:]
    labels = [str(x).lower() for x in meta.index]
    if not (
        labels[0].startswith("tumor")
        and labels[1].startswith("malignant")
        and labels[2].startswith("non-malignant")
    ):
        raise ValueError(f"Unexpected metadata rows: {list(meta.index)}")

    obs = pd.DataFrame(index=expr.columns.astype(str))
    obs["tumour"] = pd.Categorical(
        "Mel" + meta.iloc[0].astype(float).round().astype(int).astype(str).values
    )
    obs["author_malignant"] = _decode(meta.iloc[1], MALIGNANT_CODES, "malignant").to_numpy()
    obs["author_cell_type"] = _decode(meta.iloc[2], CELL_TYPE_CODES, "cell type").to_numpy()

    genes = pd.Index(expr.index.astype(str))
    n_dup = int(genes.duplicated().sum())
    adata = ad.AnnData(
        X=expr.to_numpy(dtype=np.float32).T,
        obs=obs,
        var=pd.DataFrame(index=ad.utils.make_index_unique(genes)),
    )
    adata.uns["n_duplicate_genes"] = n_dup
    adata.uns["value_scale"] = "log2(TPM/10 + 1)"
    return adata


def implied_tpm_per_cell(adata: ad.AnnData) -> np.ndarray:
    """Undo log2(TPM/10 + 1) and sum per cell.

    If the values really are log2(TPM/10 + 1), each cell should sum to close to 1,000,000
    (a little less if the authors dropped some genes before sharing).
    """
    x = np.asarray(adata.X, dtype=np.float64)
    return (10.0 * (np.power(2.0, x) - 1.0)).sum(axis=1)

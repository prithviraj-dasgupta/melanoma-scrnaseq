"""Phase 1: load the matrix and check what it actually contains before any analysis."""

from pathlib import Path

import numpy as np
import pandas as pd

from melscrna.load import implied_tpm_per_cell, read_tirosh

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "GSE72056_melanoma_single_cell_revised_v2.txt.gz"
OUT_H5AD = ROOT / "data" / "processed" / "tirosh_raw.h5ad"
OUT_TABLES = ROOT / "results" / "tables"

if __name__ == "__main__":
    adata = read_tirosh(RAW)
    x = np.asarray(adata.X)
    tpm_sums = implied_tpm_per_cell(adata)
    mt_genes = [g for g in adata.var_names if g.upper().startswith("MT-")]
    genes_per_cell = (x > 0).sum(axis=1)

    print(f"Cells: {adata.n_obs:,}   Genes: {adata.n_vars:,}")
    print(f"Duplicate gene symbols made unique: {adata.uns['n_duplicate_genes']}")
    print(f"Value range: {x.min():.2f} to {x.max():.2f}")
    print(f"All values whole numbers (would suggest raw counts): {np.allclose(x, np.round(x))}")
    print(
        "Implied TPM per cell (undoing log2(TPM/10+1)): "
        f"median {np.median(tpm_sums):,.0f}, range {tpm_sums.min():,.0f} to {tpm_sums.max():,.0f}"
    )
    print(f"Mitochondrial (MT-) genes present: {len(mt_genes)} {mt_genes[:13]}")
    print(
        f"Genes detected per cell: median {np.median(genes_per_cell):,.0f}, "
        f"min {genes_per_cell.min():,}, max {genes_per_cell.max():,}"
    )
    print(f"Tumours: {adata.obs['tumour'].nunique()}")
    print("\nAuthor malignancy labels:")
    print(adata.obs["author_malignant"].value_counts().to_string())
    print("\nAuthor cell types (non-malignant cells):")
    print(adata.obs["author_cell_type"].value_counts().to_string())

    OUT_TABLES.mkdir(parents=True, exist_ok=True)
    pd.crosstab(adata.obs["tumour"], adata.obs["author_malignant"]).to_csv(
        OUT_TABLES / "cells_per_tumour.tsv", sep="\t"
    )
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT_H5AD)
    print(f"\nSaved {OUT_H5AD.relative_to(ROOT)} and {OUT_TABLES.relative_to(ROOT)}/")

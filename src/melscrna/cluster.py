"""Highly variable genes, PCA, neighbour graph, UMAP and Leiden clustering.

The input is already log2(TPM/10 + 1), so there is no count normalisation step here
(docs/decisions.md, D3).
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def embed_and_cluster(
    adata: ad.AnnData,
    n_top_genes: int = 2000,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    resolution: float = 0.5,
    seed: int = 0,
) -> ad.AnnData:
    """Return a copy of ``adata`` with PCA, UMAP and ``obs['leiden']`` added.

    ``.X`` keeps the full log-expression matrix for marker plots; scaling is done on a
    temporary HVG-only copy.
    """
    adata = adata.copy()
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=n_top_genes)
    hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=max(n_pcs, 2), random_state=seed)
    adata.obsm["X_pca"] = hvg.obsm["X_pca"]
    adata.uns["pca"] = hvg.uns["pca"]
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(
        adata,
        resolution=resolution,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=seed,
    )
    return adata


def dominant_share(groups: pd.Series, labels: pd.Series) -> pd.Series:
    """For each group (e.g. cluster), the fraction of its cells that carry its most common label.

    1.0 means every cell in the cluster comes from one label (e.g. one patient).
    """
    table = pd.crosstab(np.asarray(groups), np.asarray(labels))
    return table.max(axis=1) / table.sum(axis=1)

"""Phase 3: HVGs, PCA, neighbours, UMAP and Leiden clustering, and what drives the clusters."""

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from melscrna.cluster import dominant_share, embed_and_cluster

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD = ROOT / "data" / "processed" / "tirosh_qc.h5ad"
OUT_H5AD = ROOT / "data" / "processed" / "tirosh_clustered.h5ad"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
N_PCS = 30
RESOLUTION = 0.5

if __name__ == "__main__":
    adata = embed_and_cluster(ad.read_h5ad(IN_H5AD), n_pcs=N_PCS, resolution=RESOLUTION)
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    vr = adata.uns["pca"]["variance_ratio"]
    print(f"HVGs: {int(adata.var['highly_variable'].sum())}")
    print(f"Variance explained by PC1-{N_PCS}: {vr[:N_PCS].sum():.1%} (PC1 alone {vr[0]:.1%})")
    print(f"Leiden clusters (resolution {RESOLUTION}): {adata.obs['leiden'].nunique()}")

    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(np.arange(1, len(vr) + 1), 100 * vr, "o-", ms=3)
    ax.axvline(N_PCS, color="grey", ls="--", lw=1)
    ax.set(
        xlabel="Principal component",
        ylabel="Variance explained (%)",
        title=f"PCA elbow plot (dashed: {N_PCS} PCs used)",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "pca_elbow.png", dpi=150)

    composition = pd.crosstab(adata.obs["leiden"], adata.obs["author_malignant"])
    composition["dominant_tumour"] = pd.crosstab(adata.obs["leiden"], adata.obs["tumour"]).idxmax(
        axis=1
    )
    composition["dominant_tumour_share"] = dominant_share(adata.obs["leiden"], adata.obs["tumour"])
    composition["n_tumours"] = (pd.crosstab(adata.obs["leiden"], adata.obs["tumour"]) > 0).sum(1)
    composition.to_csv(TABLES / "cluster_composition.tsv", sep="\t")
    print("\n" + composition.round(2).to_string())

    frac_malignant = composition["malignant"] / composition[
        ["malignant", "non-malignant", "unresolved"]
    ].sum(axis=1)
    mal = frac_malignant >= 0.8
    non = frac_malignant <= 0.2
    print(
        f"\nMostly-malignant clusters (>=80% malignant): {int(mal.sum())}; median share of cells "
        f"from one tumour = {composition.loc[mal, 'dominant_tumour_share'].median():.0%}"
    )
    print(
        f"Mostly-non-malignant clusters (<=20% malignant): {int(non.sum())}; median share "
        f"from one tumour = {composition.loc[non, 'dominant_tumour_share'].median():.0%}"
    )

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    panels = [
        ("leiden", "on data"),
        ("author_malignant", "right margin"),
        ("author_cell_type", "right margin"),
        ("tumour", "right margin"),
    ]
    for ax, (key, legend) in zip(axes.ravel(), panels, strict=True):
        sc.pl.umap(
            adata,
            color=key,
            ax=ax,
            legend_loc=legend,
            legend_fontsize=7,
            frameon=False,
            show=False,
            size=12,
        )
    fig.tight_layout()
    fig.savefig(FIGURES / "umap_overview.png", dpi=150)
    adata.write_h5ad(OUT_H5AD)
    print(f"\nSaved {OUT_H5AD.relative_to(ROOT)}")

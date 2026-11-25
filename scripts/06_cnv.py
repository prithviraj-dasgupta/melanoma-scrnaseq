"""Phase 5: infer copy-number variation with infercnvpy and call malignant cells.

The reference is half of the cells in clusters we labelled as immune or stromal (our own
annotation, not the author labels). The other half is held out to set the thresholds.
"""

from pathlib import Path

import anndata as ad
import infercnvpy as cnv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from melscrna.cnv import (
    add_gene_positions,
    call_malignant,
    cnv_correlation,
    cnv_signal,
    immune_locus_genes,
    read_gene_positions,
    split_reference,
)

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD = ROOT / "data" / "processed" / "tirosh_annotated.h5ad"
GTF = ROOT / "data" / "reference" / "gencode.v19.annotation.gtf.gz"
OUT_H5AD = ROOT / "data" / "processed" / "tirosh_cnv.h5ad"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
NORMAL_TYPES = ["T cell", "B cell", "Plasma cell", "Macrophage", "pDC", "Endothelial", "CAF"]

if __name__ == "__main__":
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(IN_H5AD)
    n_matched = add_gene_positions(adata, read_gene_positions(GTF))
    print(f"Genes with a genomic position: {n_matched:,} of {adata.n_vars:,}")

    adata.obs["cnv_role"] = split_reference(adata.obs["cell_type"].astype(str), NORMAL_TYPES)
    print(adata.obs["cnv_role"].str.split(":").str[0].value_counts().to_string())
    immune_loci = immune_locus_genes(adata.var_names)
    adata.var.loc[immune_loci, "chromosome"] = np.nan  # infercnvpy skips genes with no position
    print(f"HLA / immunoglobulin / T-cell receptor genes excluded: {int(immune_loci.sum())}")
    held_out = (adata.obs["cnv_role"] == "held_out").to_numpy()

    def run(data, per_type: bool, window: int):
        """Run infercnvpy and the malignant call for one setting."""
        key = "cnv_role" if per_type else "cnv_role_single"
        data.obs["cnv_role_single"] = data.obs["cnv_role"].str.split(":").str[0]
        cats = sorted(c for c in data.obs[key].unique() if c.startswith("reference"))
        cnv.tl.infercnv(data, reference_key=key, reference_cat=cats, window_size=window, n_jobs=1)
        sig = cnv_signal(data.obsm["X_cnv"])
        cor = cnv_correlation(data.obsm["X_cnv"], data.obs["tumour"], sig)
        return (sig, cor, *call_malignant(sig, cor, held_out, quantile=0.99))

    def score(calls_):
        obs = adata.obs
        res = (obs["author_malignant"] != "unresolved").to_numpy()
        a = (obs["author_malignant"] == "malignant").to_numpy()[res]
        o = calls_[res]
        stromal = obs["cell_type"].isin(["CAF", "Endothelial", "Macrophage"]).to_numpy()
        return {
            "agreement": round(float((a == o).mean()), 3),
            "sensitivity": round(float((a & o).sum() / a.sum()), 3),
            "specificity": round(float((~a & ~o).sum() / (~a).sum()), 3),
            "stromal_or_myeloid_called_malignant": int(calls_[stromal].sum()),
        }

    # Sensitivity analysis: how much do the results depend on two settings?
    rows = []
    for per_type in (True, False):
        for window in (100, 150):
            *_, variant_calls, _, _ = run(adata.copy(), per_type, window)
            rows.append(
                {
                    "reference": "one per cell type" if per_type else "pooled",
                    "window_genes": window,
                    **score(variant_calls),
                }
            )
    sensitivity = pd.DataFrame(rows)
    sensitivity.to_csv(TABLES / "cnv_sensitivity_analysis.tsv", sep="\t", index=False)
    print("\nSensitivity analysis (main analysis = first row):")
    print(sensitivity.to_string(index=False))

    # Main analysis: one reference per cell type, default 100-gene window.
    signal, corr, calls, s_thr, c_thr = run(adata, per_type=True, window=100)
    adata.obs["cnv_signal"] = signal
    adata.obs["cnv_correlation"] = corr
    adata.obs["cnv_call"] = pd.Categorical(np.where(calls, "malignant", "non-malignant"))
    print(
        f"\nThresholds (99th percentile of held-out normal cells): "
        f"signal {s_thr:.4f}, correlation {c_thr:.3f}"
    )
    print(
        f"Held-out normal cells called malignant: {int(calls[held_out].sum())} of "
        f"{int(held_out.sum())} ({calls[held_out].mean():.1%})"
    )

    # Comparison with the author labels (they were also partly CNV-based: see D4).
    confusion = pd.crosstab(adata.obs["author_malignant"], adata.obs["cnv_call"])
    confusion.to_csv(TABLES / "cnv_vs_author.tsv", sep="\t")
    print("\n" + confusion.to_string())
    resolved = adata.obs["author_malignant"] != "unresolved"
    a_mal = (adata.obs["author_malignant"] == "malignant")[resolved]
    o_mal = (adata.obs["cnv_call"] == "malignant")[resolved]
    sens = (a_mal & o_mal).sum() / a_mal.sum()
    spec = (~a_mal & ~o_mal).sum() / (~a_mal).sum()
    acc = (a_mal == o_mal).mean()
    print(
        f"\nVs author labels (unresolved excluded): agreement {acc:.1%}, "
        f"sensitivity {sens:.1%}, specificity {spec:.1%}"
    )

    mal_obs = adata.obs[adata.obs["author_malignant"] == "malignant"]
    per_tumour = pd.crosstab(mal_obs["tumour"], mal_obs["cnv_call"])
    per_tumour["sensitivity"] = (per_tumour["malignant"] / per_tumour.sum(axis=1)).round(3)
    per_tumour.to_csv(TABLES / "cnv_sensitivity_per_tumour.tsv", sep="\t")
    print("\nAuthor-malignant cells called malignant, per tumour:")
    print(per_tumour.to_string())

    by_type = pd.crosstab(adata.obs["cell_type"], adata.obs["cnv_call"])
    by_type.to_csv(TABLES / "cnv_by_cell_type.tsv", sep="\t")
    print("\nCNV calls by our marker-based label:")
    print(by_type.to_string())

    # Figure 1: the two CNV measures, with thresholds.
    fig, ax = plt.subplots(figsize=(6.5, 5))
    colours = {"malignant": "#c0392b", "non-malignant": "#2e86c1", "unresolved": "#7f8c8d"}
    for label, colour in colours.items():
        sel = (adata.obs["author_malignant"] == label).to_numpy()
        ax.scatter(signal[sel], corr[sel], s=4, alpha=0.5, c=colour, label=f"author: {label}")
    ax.axvline(s_thr, color="black", ls="--", lw=1)
    ax.axhline(c_thr, color="black", ls="--", lw=1)
    ax.set(
        xlabel="CNV signal (mean squared CNV)",
        ylabel="Correlation with tumour CNV profile",
        title="Malignant = beyond either dashed line",
    )
    ax.legend(frameon=False, markerscale=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "cnv_signal_vs_correlation.png", dpi=150)
    plt.close(fig)

    # Figure 2: CNV heatmap, one row group per tumour's malignant cells plus normal cells.
    adata.obs["cnv_group"] = np.where(
        adata.obs["cnv_call"] == "malignant",
        adata.obs["tumour"].astype(str),
        np.where(adata.obs["cnv_role"] == "query", "normal (query)", "normal (reference/held-out)"),
    )
    adata.obs["cnv_group"] = adata.obs["cnv_group"].astype("category")
    cnv.pl.chromosome_heatmap(adata, groupby="cnv_group", show=False, figsize=(12, 9))
    plt.savefig(FIGURES / "cnv_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close("all")

    # Figure 3: our CNV call next to the author label on the Phase 3 UMAP.
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, key, title in [
        (axes[0], "cnv_call", "Our CNV-based call"),
        (axes[1], "author_malignant", "Author label"),
    ]:
        sc.pl.umap(
            adata,
            color=key,
            ax=ax,
            frameon=False,
            show=False,
            size=10,
            title=title,
            palette={"malignant": "#c0392b", "non-malignant": "#2e86c1", "unresolved": "#95a5a6"},
        )
    fig.tight_layout()
    fig.savefig(FIGURES / "umap_cnv_call.png", dpi=150)
    plt.close(fig)

    adata.write_h5ad(OUT_H5AD)
    print(f"\nSaved {OUT_H5AD.relative_to(ROOT)}")

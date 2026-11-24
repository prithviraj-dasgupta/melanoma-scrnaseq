"""Phase 2: quality control. Reports how many cells and genes each step removes, and why."""

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from melscrna.qc import MELANOCYTE_MARKERS, add_qc_metrics, mad_bounds, run_qc

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD = ROOT / "data" / "processed" / "tirosh_raw.h5ad"
OUT_H5AD = ROOT / "data" / "processed" / "tirosh_qc.h5ad"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
N_MADS = 5.0

if __name__ == "__main__":
    raw = ad.read_h5ad(IN_H5AD)
    result = run_qc(raw, n_mads=N_MADS, conflict_threshold=7.0, min_cells=3)
    summary = result.summary()
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    summary.to_csv(TABLES / "qc_summary.tsv", sep="\t", index=False)
    with pd.option_context("display.max_colwidth", 90, "display.width", 200):
        print(summary.to_string(index=False))

    # What a stricter, more "standard" 3-MAD cut-off would have done, and to whom.
    add_qc_metrics(raw)
    lo3, hi3 = mad_bounds(raw.obs["n_genes"], n_mads=3)
    strict_out = (raw.obs["n_genes"] < lo3) | (raw.obs["n_genes"] > hi3)
    print(
        f"\nFor comparison, a 3-MAD cut-off ({lo3:,.0f} to {hi3:,.0f}) would remove "
        f"{int(strict_out.sum())} cells:"
    )
    print(raw.obs.loc[strict_out, "author_malignant"].value_counts().to_string())
    print("\nMedian genes detected by author label:")
    print(raw.obs.groupby("author_malignant", observed=True)["n_genes"].median().to_string())

    # Low-level melanocyte transcripts in immune cells (reported, not filtered).
    x = raw[:, list(MELANOCYTE_MARKERS)].X
    mel = np.asarray(x).max(axis=1)
    t_cells = (raw.obs["author_cell_type"] == "T cell").to_numpy()
    print(
        f"\nT cells with any PMEL/MLANA/TYR > 4 (~150 TPM): "
        f"{int((mel[t_cells] > 4).sum())} of {int(t_cells.sum())}"
    )

    low, high = result.adata.uns["qc_n_genes_bounds"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for label, colour in [
        ("malignant", "#c0392b"),
        ("non-malignant", "#2e86c1"),
        ("unresolved", "#7f8c8d"),
    ]:
        sel = raw.obs["author_malignant"] == label
        axes[0].hist(raw.obs.loc[sel, "n_genes"], bins=60, alpha=0.6, label=label, color=colour)
        axes[1].hist(raw.obs.loc[sel, "pct_ribo"], bins=60, alpha=0.6, label=label, color=colour)
    for v in (low, high):
        axes[0].axvline(v, color="black", ls="--", lw=1)
    for v in (lo3, hi3):
        axes[0].axvline(v, color="grey", ls=":", lw=1)
    axes[0].set(
        xlabel="Genes detected per cell",
        ylabel="Cells",
        title="Genes detected (dashed: 5-MAD bounds used; dotted: 3-MAD)",
    )
    axes[1].set(
        xlabel="Ribosomal protein genes (% of TPM)",
        ylabel="Cells",
        title="Ribosomal share (reported, not filtered)",
    )
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "qc_distributions.png", dpi=150)

    result.adata.write_h5ad(OUT_H5AD)
    print(
        f"\nKept {result.adata.n_obs:,} cells x {result.adata.n_vars:,} genes -> "
        f"{OUT_H5AD.relative_to(ROOT)}"
    )

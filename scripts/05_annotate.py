"""Phase 4: label clusters from marker genes, back each label with data-driven top genes,
and compare with the author cell-type labels."""

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from melscrna.annotate import MARKERS, TO_AUTHOR, label_from_scores, marker_scores

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD = ROOT / "data" / "processed" / "tirosh_clustered.h5ad"
OUT_H5AD = ROOT / "data" / "processed" / "tirosh_annotated.h5ad"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"

if __name__ == "__main__":
    adata = ad.read_h5ad(IN_H5AD)
    scores = marker_scores(adata, "leiden")
    labels = label_from_scores(scores)

    # Data-driven check: top genes for each cluster vs the rest (Wilcoxon rank-sum).
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    labels["top_genes"] = [
        ", ".join(sc.get.rank_genes_groups_df(adata, str(c)).names[:8]) for c in labels.index
    ]
    labels["n_cells"] = adata.obs["leiden"].value_counts().reindex(labels.index.astype(str)).values
    labels.index = labels.index.astype(str)
    TABLES.mkdir(parents=True, exist_ok=True)
    labels.to_csv(TABLES / "cluster_annotation.tsv", sep="\t")
    scores.round(2).to_csv(TABLES / "cluster_marker_scores.tsv", sep="\t")
    with pd.option_context("display.width", 250, "display.max_colwidth", 70):
        print(
            labels[
                ["n_cells", "label", "top_score", "second", "second_score", "top_genes"]
            ].to_string()
        )

    adata.obs["cell_type"] = adata.obs["leiden"].astype(str).map(labels["label"]).astype("category")

    # Agreement with the author cell types, for cells that have one.
    has_type = adata.obs["author_cell_type"] != "none"
    ours = adata.obs.loc[has_type, "cell_type"].astype(str).map(TO_AUTHOR).fillna("other")
    theirs = adata.obs.loc[has_type, "author_cell_type"].astype(str)
    agree = ours == theirs
    print(
        f"\nAgreement with author cell types: {int(agree.sum()):,} of {int(has_type.sum()):,} "
        f"cells ({agree.mean():.1%})"
    )
    per_type = pd.DataFrame({"author": theirs, "agree": agree}).groupby("author")["agree"]
    print(per_type.agg(["sum", "count", "mean"]).round(3).to_string())
    pd.crosstab(theirs, adata.obs.loc[has_type, "cell_type"]).to_csv(
        TABLES / "annotation_vs_author.tsv", sep="\t"
    )

    mal = adata.obs["author_malignant"] == "malignant"
    melano = adata.obs["cell_type"].astype(str).str.contains("Melanocytic")
    print(
        f"\nAuthor-malignant cells in a Melanocytic (or mixed Melanocytic) cluster: "
        f"{int((mal & melano).sum())} of {int(mal.sum())}"
    )
    print("Where the other author-malignant cells sit:")
    print(adata.obs.loc[mal & ~melano, "cell_type"].value_counts().loc[lambda s: s > 0].to_string())

    FIGURES.mkdir(parents=True, exist_ok=True)
    marker_genes = {k: [g for g in v if g in adata.var_names] for k, v in MARKERS.items()}
    sc.pl.dotplot(adata, marker_genes, groupby="leiden", standard_scale="var", show=False)
    plt.savefig(FIGURES / "marker_dotplot.png", dpi=150, bbox_inches="tight")
    plt.close("all")
    fig, ax = plt.subplots(figsize=(9, 7))
    sc.pl.umap(
        adata,
        color="cell_type",
        ax=ax,
        frameon=False,
        show=False,
        size=12,
        title="Marker-based annotation (per cluster)",
    )
    fig.savefig(FIGURES / "umap_annotation.png", dpi=150, bbox_inches="tight")

    adata.write_h5ad(OUT_H5AD)
    print(f"\nSaved {OUT_H5AD.relative_to(ROOT)}")

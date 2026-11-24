"""Label clusters from known marker genes.

Each cell type has a short list of canonical markers. For every cluster we take the mean
expression of each list, and the cluster gets the highest-scoring cell type, unless the
evidence is weak (``Unclear``) or two cell types score highly (``Mixed``).

"Melanocytic" means melanocyte-lineage genes are on. It does **not** mean malignant: that
call is made from inferred copy-number changes (Phase 5).
"""

from __future__ import annotations

import warnings

import anndata as ad
import numpy as np
import pandas as pd

MARKERS: dict[str, list[str]] = {
    "T cell": ["CD3D", "CD3E", "CD2"],
    "NK cell": ["NCAM1", "KLRF1", "KLRC1"],
    "B cell": ["CD19", "MS4A1", "CD79A"],
    # MZB1 and IGJ (JCHAIN) are also high in pDCs here, so plasma cells use CD138, BLIMP1, BCMA
    "Plasma cell": ["SDC1", "PRDM1", "TNFRSF17"],
    "Macrophage": ["CD14", "CD68", "CD163", "AIF1"],
    "pDC": ["LILRA4", "IL3RA", "CLEC4C", "PLD4"],  # plasmacytoid dendritic cells
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "CAF": ["COL1A1", "COL1A2", "COL3A1", "DCN"],  # cancer-associated fibroblasts
    "Melanocytic": ["PMEL", "MLANA", "TYR", "MITF", "SOX10"],
}

# Pairs that share genes for biological reasons, so scoring high on both is expected, not
# a sign of a mixed cluster: plasma cells develop from B cells and keep CD79A; NK and
# cytotoxic T cells share killing machinery.
RELATED = {frozenset({"Plasma cell", "B cell"}), frozenset({"NK cell", "T cell"})}

# Our labels mapped onto the author categories, for the agreement check.
TO_AUTHOR = {
    "T cell": "T cell",
    "NK cell": "NK cell",
    "B cell": "B cell",
    "Plasma cell": "B cell",
    "Macrophage": "Macrophage",
    "Endothelial": "Endothelial",
    "CAF": "CAF",
}


def marker_scores(
    adata: ad.AnnData, groupby: str, markers: dict[str, list[str]] = MARKERS
) -> pd.DataFrame:
    """Mean expression of each marker list, averaged over the cells of each group.

    Returns a groups x cell types table. Marker genes absent from the data are skipped with a
    warning; a cell type with no markers present raises ``KeyError``.
    """
    per_cell = {}
    for cell_type, genes in markers.items():
        present = [g for g in genes if g in adata.var_names]
        missing = sorted(set(genes) - set(present))
        if not present:
            raise KeyError(f"No markers for {cell_type} are in the data: {genes}")
        if missing:
            warnings.warn(f"{cell_type}: markers not in data: {missing}", stacklevel=2)
        per_cell[cell_type] = np.asarray(adata[:, present].X).mean(axis=1)
    scores = pd.DataFrame(per_cell, index=adata.obs_names)
    return scores.groupby(adata.obs[groupby].to_numpy()).mean()


def label_from_scores(
    scores: pd.DataFrame,
    min_score: float = 2.5,
    mixed_score: float = 2.0,
    mixed_ratio: float = 0.5,
) -> pd.DataFrame:
    """Turn a groups x cell types score table into one label per group.

    - top score below ``min_score`` -> ``Unclear`` (no marker list is convincingly on)
    - second-best score at least ``mixed_score`` and at least ``mixed_ratio`` x the top score
      -> ``Mixed: A / B`` (two lineages clearly on in the same cluster), unless the two are
      in ``RELATED``
    - otherwise the top cell type
    """
    rows = []
    for group, row in scores.iterrows():
        ranked = row.sort_values(ascending=False)
        top, second = ranked.index[0], ranked.index[1]
        if ranked.iloc[0] < min_score:
            label = "Unclear"
        elif ranked.iloc[1] >= max(mixed_score, mixed_ratio * ranked.iloc[0]) and (
            frozenset({top, second}) not in RELATED
        ):
            label = f"Mixed: {top} / {second}"
        else:
            label = top
        rows.append(
            {
                "group": group,
                "label": label,
                "top": top,
                "top_score": round(float(ranked.iloc[0]), 2),
                "second": second,
                "second_score": round(float(ranked.iloc[1]), 2),
            }
        )
    return pd.DataFrame(rows).set_index("group")

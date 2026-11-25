"""Infer copy-number variation (CNV) from expression and call malignant cells.

Idea: a chromosome arm that is gained makes many neighbouring genes slightly more expressed;
a loss makes them slightly less. Averaging expression over windows of ~100 neighbouring
genes, relative to normal reference cells, turns that into a CNV profile per cell. Normal
cells have flat profiles; tumour cells carry their tumour's gains and losses.

Malignant calls follow the two measures used by Tirosh et al. (2016):

- ``cnv_signal``: mean squared CNV value (how far the profile is from flat)
- ``cnv_correlation``: correlation with the average profile of the strongest CNV cells in
  the same tumour (does the cell share its tumour's specific changes)

Thresholds come from held-out normal cells that were *not* used as the reference, so the
cut-off is not biased by the reference being flat by construction.
"""

from __future__ import annotations

import gzip
from collections.abc import Sequence
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

MAIN_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def read_gene_positions(gtf_path: str | Path) -> pd.DataFrame:
    """Read gene coordinates from a GENCODE GTF: one row per gene name.

    Only ``gene`` records on chr1-22 and chrX are kept (chrY is dropped because it would
    track patient sex, not tumour copy number). If a name appears more than once, the first
    record is kept.
    """
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    rows = []
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene" or fields[0] not in MAIN_CHROMOSOMES:
                continue
            attrs = fields[8]
            key = 'gene_name "'
            i = attrs.find(key)
            if i < 0:
                continue
            name = attrs[i + len(key) : attrs.index('"', i + len(key))]
            rows.append((name, fields[0], int(fields[3]), int(fields[4])))
    table = pd.DataFrame(rows, columns=["gene", "chromosome", "start", "end"])
    return table.drop_duplicates("gene").set_index("gene")


def add_gene_positions(adata: ad.AnnData, positions: pd.DataFrame) -> int:
    """Add chromosome/start/end to ``adata.var`` in place; return how many genes matched."""
    matched = positions.reindex(adata.var_names)
    adata.var["chromosome"] = matched["chromosome"].to_numpy()
    adata.var["start"] = matched["start"].to_numpy()
    adata.var["end"] = matched["end"].to_numpy()
    return int(matched["chromosome"].notna().sum())


# Gene families whose expression is specific to immune cells and clustered in the genome
# (MHC/HLA on chr6p21, immunoglobulin and T-cell receptor loci). Against a reference of
# immune cells they look like deletions in every non-immune cell, so they are excluded.
IMMUNE_LOCI_PATTERN = (
    r"^(HLA-|IGH[VDJ]\d|IGH(G\d|A\d|M|E|D)$|IG[KL][VJC]|TR[ABDG](V\d|D\d|J\d|C\d?$))"
)


def immune_locus_genes(var_names: pd.Index) -> np.ndarray:
    """Boolean mask of HLA, immunoglobulin and T-cell receptor genes."""
    return np.asarray(pd.Index(var_names).str.match(IMMUNE_LOCI_PATTERN), dtype=bool)


def split_reference(
    labels: pd.Series, normal_types: Sequence[str], fraction: float = 0.5, seed: int = 0
) -> pd.Series:
    """Split normal cells into ``reference`` and ``held_out`` halves, within each cell type.

    Everything else is ``query``. The reference defines "flat"; the held-out normal cells
    give an honest estimate of how much CNV signal a normal cell shows. Reference cells are
    labelled ``reference: <cell type>`` so infercnvpy can use each type as its own baseline.
    """
    rng = np.random.default_rng(seed)
    role = pd.Series("query", index=labels.index, dtype=object)
    for cell_type in normal_types:
        idx = labels.index[labels == cell_type].to_numpy()
        if len(idx) == 0:
            continue
        chosen = rng.permutation(idx)[: max(1, round(fraction * len(idx)))]
        role.loc[idx] = "held_out"
        role.loc[chosen] = f"reference: {cell_type}"
    return role


def cnv_signal(cnv_matrix: np.ndarray) -> np.ndarray:
    """Mean squared CNV value per cell."""
    m = np.asarray(cnv_matrix.todense() if hasattr(cnv_matrix, "todense") else cnv_matrix)
    return (m**2).mean(axis=1)


def cnv_correlation(
    cnv_matrix: np.ndarray, groups: Sequence, signal: np.ndarray, top_fraction: float = 0.1
) -> np.ndarray:
    """Correlation of each cell's CNV profile with its group's (tumour's) consensus profile.

    The consensus is the mean profile of the ``top_fraction`` of cells in that group with the
    highest ``signal``: the cells most likely to be malignant, chosen without using any
    labels.
    """
    m = np.asarray(cnv_matrix.todense() if hasattr(cnv_matrix, "todense") else cnv_matrix)
    groups = np.asarray(groups)
    out = np.full(m.shape[0], np.nan)
    for g in np.unique(groups):
        idx = np.flatnonzero(groups == g)
        n_top = max(1, int(round(top_fraction * len(idx))))
        top = idx[np.argsort(signal[idx])[::-1][:n_top]]
        consensus = m[top].mean(axis=0)
        centred = m[idx] - m[idx].mean(axis=1, keepdims=True)
        c = consensus - consensus.mean()
        denom = np.linalg.norm(centred, axis=1) * np.linalg.norm(c)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[idx] = centred @ c / denom
    return out


def call_malignant(
    signal: np.ndarray,
    correlation: np.ndarray,
    held_out: np.ndarray,
    quantile: float = 0.99,
) -> tuple[np.ndarray, float, float]:
    """Malignant if signal OR correlation exceeds the ``quantile`` of held-out normal cells.

    Returns (calls, signal_threshold, correlation_threshold).
    """
    held_out = np.asarray(held_out, dtype=bool)
    if not held_out.any():
        raise ValueError("Need at least one held-out normal cell to set thresholds")
    s_thr = float(np.quantile(signal[held_out], quantile))
    c_thr = float(np.nanquantile(correlation[held_out], quantile))
    calls = (signal > s_thr) | (np.nan_to_num(correlation, nan=-1.0) > c_thr)
    return calls, s_thr, c_thr

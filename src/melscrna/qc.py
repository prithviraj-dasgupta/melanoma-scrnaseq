"""Quality control for log2(TPM/10 + 1) Smart-seq2 data.

Standard droplet (10x) QC uses raw counts, mitochondrial share and scrublet. None of those fit
this dataset (see docs/decisions.md, D3 and D6), so QC here is:

1. remove cells with outlying numbers of detected genes (median +/- n MADs on a log scale)
2. remove cells that strongly express both an immune and a melanocyte marker (a biological
   doublet check that works on any normalised data)
3. remove genes detected in fewer than ``min_cells`` cells

Each step records how many cells or genes it removed.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import anndata as ad
import numpy as np
import pandas as pd

IMMUNE_MARKERS = ("PTPRC",)  # CD45, on every immune cell
MELANOCYTE_MARKERS = ("PMEL", "MLANA", "TYR")  # melanocyte lineage, high in melanoma cells


def to_tpm(x: np.ndarray) -> np.ndarray:
    """Undo log2(TPM/10 + 1)."""
    return 10.0 * (np.power(2.0, np.asarray(x, dtype=np.float64)) - 1.0)


def add_qc_metrics(adata: ad.AnnData) -> None:
    """Add per-cell QC columns to ``adata.obs`` in place.

    - ``n_genes``: genes with any expression
    - ``pct_ribo``: share of TPM from ribosomal protein genes (RPS*/RPL*)
    """
    x = np.asarray(adata.X)
    adata.obs["n_genes"] = (x > 0).sum(axis=1)
    ribo = adata.var_names.str.match(r"^RP[SL]\d")
    tpm = to_tpm(x)
    adata.obs["pct_ribo"] = 100.0 * tpm[:, ribo].sum(axis=1) / tpm.sum(axis=1)


def mad_bounds(values: Sequence[float], n_mads: float, log: bool = True) -> tuple[float, float]:
    """Return (low, high) = median -/+ n_mads x median absolute deviation.

    MAD is used rather than the standard deviation because it is not pulled around by the
    outliers we are trying to find. With ``log=True`` the bounds are computed on log1p values
    and converted back.
    """
    v = np.asarray(values, dtype=float)
    if log:
        v = np.log1p(v)
    med = np.median(v)
    mad = np.median(np.abs(v - med))
    low, high = med - n_mads * mad, med + n_mads * mad
    if log:
        return float(np.expm1(low)), float(np.expm1(high))
    return float(low), float(high)


def _max_expression(adata: ad.AnnData, genes: Sequence[str]) -> np.ndarray:
    present = [g for g in genes if g in adata.var_names]
    if not present:
        raise KeyError(f"None of the marker genes are in the data: {list(genes)}")
    return np.asarray(adata[:, present].X).max(axis=1)


def flag_lineage_conflict(
    adata: ad.AnnData,
    threshold: float,
    immune: Sequence[str] = IMMUNE_MARKERS,
    melanocyte: Sequence[str] = MELANOCYTE_MARKERS,
) -> np.ndarray:
    """True for cells expressing both an immune and a melanocyte marker above ``threshold``.

    One cell cannot be both a CD45+ immune cell and a melanocyte, so strong co-expression
    suggests two cells captured together (a doublet).
    """
    return (_max_expression(adata, immune) > threshold) & (
        _max_expression(adata, melanocyte) > threshold
    )


@dataclass
class QCResult:
    adata: ad.AnnData
    steps: list[dict] = field(default_factory=list)

    def summary(self) -> pd.DataFrame:
        return pd.DataFrame(self.steps)


def run_qc(
    adata: ad.AnnData,
    n_mads: float = 5.0,
    conflict_threshold: float = 7.0,
    min_cells: int = 3,
) -> QCResult:
    """Run the three QC steps and record what each removed. Does not modify the input."""
    adata = adata.copy()
    add_qc_metrics(adata)
    steps: list[dict] = []

    def record(step: str, unit: str, removed: int, remaining: int, reason: str) -> None:
        steps.append(
            {
                "step": step,
                "unit": unit,
                "removed": removed,
                "remaining": remaining,
                "reason": reason,
            }
        )

    record("start", "cells", 0, adata.n_obs, "cells in the deposited file")

    low, high = mad_bounds(adata.obs["n_genes"], n_mads=n_mads)
    keep = ((adata.obs["n_genes"] >= low) & (adata.obs["n_genes"] <= high)).to_numpy()
    adata.uns["qc_n_genes_bounds"] = [low, high]
    record(
        "genes detected outlier",
        "cells",
        int((~keep).sum()),
        int(keep.sum()),
        f"outside {low:,.0f} to {high:,.0f} detected genes (median +/- {n_mads:g} MADs, log scale)",
    )
    adata = adata[keep].copy()

    conflict = flag_lineage_conflict(adata, threshold=conflict_threshold)
    record(
        "immune + melanocyte co-expression",
        "cells",
        int(conflict.sum()),
        int((~conflict).sum()),
        f"PTPRC and PMEL/MLANA/TYR both > {conflict_threshold:g} log2(TPM/10+1), likely doublets",
    )
    adata = adata[~conflict].copy()

    n_cells_per_gene = (np.asarray(adata.X) > 0).sum(axis=0)
    keep_genes = n_cells_per_gene >= min_cells
    steps.append(
        {
            "step": "rare genes",
            "unit": "genes",
            "removed": int((~keep_genes).sum()),
            "remaining": int(keep_genes.sum()),
            "reason": f"detected in fewer than {min_cells} cells (no information for clustering)",
        }
    )
    adata = adata[:, keep_genes].copy()
    return QCResult(adata=adata, steps=steps)

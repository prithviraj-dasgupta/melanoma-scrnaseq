"""MITF-high vs AXL-high malignant cell states.

Melanoma cells can sit in a differentiated, MITF-high state (pigmentation genes on) or an
undifferentiated, AXL-high state (invasive genes on), and AXL-high states have been linked to
resistance to BRAF/MEK inhibitors. Single genes are noisy in single-cell data, so each state
is measured as a *programme*: the genes that most track MITF or AXL.

To avoid testing on the same cells used to choose the genes, programmes are learned on one
random half of the malignant cells and scored on the other half.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


def split_half(groups: pd.Series, seed: int = 0) -> pd.Series:
    """True for a random half of the cells within each group (e.g. each tumour)."""
    rng = np.random.default_rng(seed)
    mask = pd.Series(False, index=groups.index)
    for _, idx in groups.groupby(groups, observed=True).groups.items():
        shuffled = rng.permutation(np.asarray(idx))
        mask.loc[shuffled[: len(shuffled) // 2]] = True
    return mask


def centre_within(data: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """Subtract each group's mean, so correlations reflect variation *within* tumours.

    Without this, two genes can look correlated only because some tumours are high in both.
    """
    return data - data.groupby(np.asarray(groups)).transform("mean")


def derive_program(
    centred: pd.DataFrame, anchor: str, n_genes: int = 50, exclude: Sequence[str] = ()
) -> list[str]:
    """The ``n_genes`` genes most positively correlated with ``anchor`` (anchor excluded)."""
    corr = centred.corrwith(centred[anchor]).drop(labels=[anchor, *exclude], errors="ignore")
    return list(corr.dropna().nlargest(n_genes).index)


def score_program(data: pd.DataFrame, genes: Sequence[str], reference: pd.DataFrame) -> pd.Series:
    """Mean z-score of the programme genes, using ``reference`` for each gene's mean and SD."""
    genes = list(genes)
    mu = reference[genes].mean()
    sd = reference[genes].std().replace(0, np.nan)
    return ((data[genes] - mu) / sd).mean(axis=1)

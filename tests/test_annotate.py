import anndata as ad
import numpy as np
import pandas as pd
import pytest

from melscrna.annotate import label_from_scores, marker_scores


def test_label_rules():
    scores = pd.DataFrame(
        {"T cell": [6.0, 1.0, 4.0, 8.0], "B cell": [0.5, 1.5, 3.0, 2.5]},
        index=["clear", "weak", "mixed", "related"],
    )
    labels = label_from_scores(scores, min_score=2.5, mixed_score=2.0, mixed_ratio=0.5)["label"]
    assert labels["clear"] == "T cell"
    assert labels["weak"] == "Unclear"
    assert labels["mixed"] == "Mixed: T cell / B cell"
    assert labels["related"] == "T cell"  # second score below half of the top


def test_related_lineages_are_not_called_mixed():
    scores = pd.DataFrame({"Plasma cell": [4.0], "B cell": [3.0]}, index=["plasma"])
    assert label_from_scores(scores).loc["plasma", "label"] == "Plasma cell"


def test_marker_scores_average_genes_then_cells():
    x = np.array([[4.0, 2.0, 0.0], [6.0, 4.0, 0.0], [0.0, 0.0, 8.0]], dtype=np.float32)
    adata = ad.AnnData(
        X=x,
        obs=pd.DataFrame({"cluster": ["a", "a", "b"]}, index=["c1", "c2", "c3"]),
        var=pd.DataFrame(index=["CD3D", "CD3E", "MS4A1"]),
    )
    markers = {"T cell": ["CD3D", "CD3E"], "B cell": ["MS4A1", "NOT_A_GENE"]}
    with pytest.warns(UserWarning, match="NOT_A_GENE"):
        scores = marker_scores(adata, "cluster", markers)
    assert scores.loc["a", "T cell"] == pytest.approx(4.0)  # mean of 3 and 5
    assert scores.loc["b", "B cell"] == pytest.approx(8.0)


def test_marker_scores_raises_when_no_marker_present():
    adata = ad.AnnData(
        X=np.ones((1, 1), dtype=np.float32),
        obs=pd.DataFrame({"cluster": ["a"]}, index=["c1"]),
        var=pd.DataFrame(index=["CD3D"]),
    )
    with pytest.raises(KeyError):
        marker_scores(adata, "cluster", {"B cell": ["MS4A1"]})

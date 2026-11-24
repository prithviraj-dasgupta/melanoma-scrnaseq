import anndata as ad
import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import adjusted_rand_score

from melscrna.cluster import dominant_share, embed_and_cluster


def test_dominant_share():
    groups = pd.Series(["a", "a", "a", "a", "b", "b"])
    labels = pd.Series(["p1", "p1", "p1", "p2", "p3", "p3"])
    share = dominant_share(groups, labels)
    assert share["a"] == pytest.approx(0.75)
    assert share["b"] == pytest.approx(1.0)


@pytest.mark.slow
def test_embed_and_cluster_recovers_two_obvious_groups():
    rng = np.random.default_rng(0)
    n, g = 60, 200
    x = rng.gamma(1.0, 0.5, size=(n, g))
    x[:30, :50] += 6.0  # group A expresses genes 0-49
    x[30:, 50:100] += 6.0  # group B expresses genes 50-99
    adata = ad.AnnData(
        X=x.astype(np.float32),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(g)]),
    )
    out = embed_and_cluster(adata, n_top_genes=150, n_pcs=10, n_neighbors=10, resolution=0.3)
    truth = ["A"] * 30 + ["B"] * 30
    assert adjusted_rand_score(truth, out.obs["leiden"]) == pytest.approx(1.0)
    assert out.obsm["X_umap"].shape == (n, 2)
    assert out.shape == adata.shape  # full gene matrix kept

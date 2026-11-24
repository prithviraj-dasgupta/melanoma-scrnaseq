import anndata as ad
import numpy as np
import pandas as pd
import pytest

from melscrna.qc import add_qc_metrics, flag_lineage_conflict, mad_bounds, run_qc


def make_adata(x, genes):
    x = np.asarray(x, dtype=np.float32)
    return ad.AnnData(
        X=x,
        obs=pd.DataFrame(index=[f"c{i}" for i in range(x.shape[0])]),
        var=pd.DataFrame(index=genes),
    )


def test_mad_bounds_without_log():
    low, high = mad_bounds([1, 2, 3, 4, 100], n_mads=2, log=False)
    # median 3, absolute deviations [2,1,0,1,97] -> MAD 1
    assert (low, high) == (1.0, 5.0)


def test_mad_bounds_log_round_trips_the_median():
    low, high = mad_bounds([100, 100, 100], n_mads=3)
    assert low == pytest.approx(100) and high == pytest.approx(100)


def test_add_qc_metrics_counts_genes_and_ribosomal_share():
    one_k = np.log2(1000 / 10 + 1)  # 1,000 TPM
    adata = make_adata([[one_k, one_k, 0.0], [one_k, 0.0, 0.0]], ["RPL3", "CD3E", "MITF"])
    add_qc_metrics(adata)
    assert list(adata.obs["n_genes"]) == [2, 1]
    assert adata.obs["pct_ribo"].tolist() == pytest.approx([50.0, 100.0])


def test_lineage_conflict_needs_both_markers_high():
    adata = make_adata(
        [[9, 9, 0], [9, 0, 0], [0, 9, 0], [9, 3, 0]],
        ["PTPRC", "PMEL", "CD3E"],
    )
    assert flag_lineage_conflict(adata, threshold=7).tolist() == [True, False, False, False]


def test_lineage_conflict_missing_markers_raise():
    adata = make_adata([[1.0]], ["CD3E"])
    with pytest.raises(KeyError):
        flag_lineage_conflict(adata, threshold=7)


def test_run_qc_removes_and_reports_each_step():
    rng = np.random.default_rng(0)
    n_cells, n_genes = 40, 60
    x = (rng.random((n_cells, n_genes)) < 0.5) * 5.0
    genes = ["PTPRC", "PMEL"] + [f"G{i}" for i in range(n_genes - 3)] + ["RARE"]
    x[:, 0] = 0.0
    x[:, 1] = 0.0
    x[:, -1] = 0.0
    x[0, -1] = 5.0  # RARE seen in one cell only
    x[1, :] = 0.0
    x[1, 2] = 5.0  # cell 1: one gene detected -> low outlier
    x[2, 0] = x[2, 1] = 9.0  # cell 2: immune + melanocyte -> conflict
    adata = make_adata(x, genes)

    result = run_qc(adata, n_mads=10, conflict_threshold=7, min_cells=3)
    summary = result.summary().set_index("step")

    assert summary.loc["genes detected outlier", "removed"] == 1
    assert summary.loc["immune + melanocyte co-expression", "removed"] == 1
    assert result.adata.n_obs == n_cells - 2
    assert "RARE" not in result.adata.var_names
    assert "c1" not in result.adata.obs_names and "c2" not in result.adata.obs_names
    assert adata.n_obs == n_cells  # input untouched

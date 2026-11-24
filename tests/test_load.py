import gzip

import numpy as np
import pytest

from melscrna.load import implied_tpm_per_cell, read_tirosh


def test_shape_is_cells_by_genes(tiny_tirosh):
    adata = read_tirosh(tiny_tirosh)
    assert adata.shape == (4, 5)
    assert list(adata.obs_names) == ["c1", "c2", "c3", "c4"]


def test_author_labels_are_decoded(tiny_tirosh):
    obs = read_tirosh(tiny_tirosh).obs
    assert list(obs["tumour"]) == ["Mel79", "Mel79", "Mel80", "Mel80"]
    assert list(obs["author_malignant"]) == [
        "malignant",
        "non-malignant",
        "non-malignant",
        "unresolved",
    ]
    assert list(obs["author_cell_type"]) == ["none", "T cell", "Macrophage", "none"]


def test_values_are_not_changed(tiny_tirosh):
    adata = read_tirosh(tiny_tirosh)
    assert adata[0, "MITF"].X.item() == pytest.approx(8.1)
    assert adata.uns["value_scale"] == "log2(TPM/10 + 1)"


def test_duplicate_genes_made_unique_and_counted(tiny_tirosh):
    adata = read_tirosh(tiny_tirosh)
    assert adata.var_names.is_unique
    assert "MITF-1" in adata.var_names
    assert adata.uns["n_duplicate_genes"] == 1


def test_unknown_label_code_raises(tmp_path):
    path = tmp_path / "bad.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(
            "Cell\tc1\n"
            "tumor\t1\n"
            "malignant(1=no,2=yes,0=unresolved)\t7\n"
            "non-malignant cell type\t0\n"
            "MITF\t1\n"
        )
    with pytest.raises(ValueError, match="Unexpected codes"):
        read_tirosh(path)


def test_missing_metadata_rows_raise(tmp_path):
    path = tmp_path / "nometa.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("Cell\tc1\nMITF\t1\nCD3E\t2\nCD68\t0\nSOX10\t4\n")
    with pytest.raises(ValueError, match="Unexpected metadata rows"):
        read_tirosh(path)


def test_implied_tpm_inverts_the_log_transform(tiny_tirosh):
    adata = read_tirosh(tiny_tirosh)
    adata.X[:] = np.log2(1000.0 / 10.0 + 1.0)  # every gene at 1,000 TPM
    assert implied_tpm_per_cell(adata) == pytest.approx(np.full(4, 5000.0))

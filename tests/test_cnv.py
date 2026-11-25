import gzip

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from melscrna.cnv import (
    add_gene_positions,
    call_malignant,
    cnv_correlation,
    cnv_signal,
    immune_locus_genes,
    read_gene_positions,
    split_reference,
)

GTF = (
    "##comment\n"
    'chr1\tHAVANA\tgene\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "AAA";\n'
    'chr1\tHAVANA\ttranscript\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "AAA";\n'
    'chr2\tHAVANA\tgene\t50\t80\t.\t-\t.\tgene_id "G2"; gene_name "BBB";\n'
    'chrY\tHAVANA\tgene\t10\t20\t.\t+\t.\tgene_id "G3"; gene_name "YYY";\n'
    'chr3\tHAVANA\tgene\t5\t9\t.\t+\t.\tgene_id "G4"; gene_name "AAA";\n'
)


def test_read_gene_positions_keeps_genes_on_main_chromosomes(tmp_path):
    path = tmp_path / "tiny.gtf.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(GTF)
    pos = read_gene_positions(path)
    assert list(pos.index) == ["AAA", "BBB"]  # chrY dropped, duplicate AAA keeps first
    assert pos.loc["AAA"].tolist() == ["chr1", 100, 200]


def test_add_gene_positions_counts_matches():
    adata = ad.AnnData(
        np.zeros((1, 3), dtype=np.float32), var=pd.DataFrame(index=["AAA", "BBB", "ZZZ"])
    )
    pos = pd.DataFrame({"chromosome": ["chr1"], "start": [1], "end": [2]}, index=["AAA"])
    assert add_gene_positions(adata, pos) == 1
    assert pd.isna(adata.var.loc["ZZZ", "chromosome"])


def test_immune_locus_genes():
    names = pd.Index(
        [
            "HLA-A",
            "IGHG1",
            "IGHM",
            "IGKC",
            "IGHV1-2",
            "TRBC1",
            "TRAF1",
            "TRADD",
            "IGHMBP2",
            "IGF1",
            "CD3E",
        ]
    )
    expected = [True, True, True, True, True, True, False, False, False, False, False]
    assert immune_locus_genes(names).tolist() == expected


def test_split_reference_halves_each_normal_type():
    labels = pd.Series(["T cell"] * 10 + ["CAF"] * 4 + ["Melanocytic"] * 5)
    role = split_reference(labels, ["T cell", "CAF"], fraction=0.5, seed=1)
    assert (role == "reference: T cell").sum() == 5
    assert (role == "reference: CAF").sum() == 2
    assert (role == "held_out").sum() == 7
    assert (role[labels == "Melanocytic"] == "query").all()


def test_cnv_signal_and_correlation():
    profile = np.array([1.0, 1.0, -1.0, -1.0, 0.0, 0.0])
    m = np.vstack([profile, profile * 0.5, np.zeros(6) + 0.01 * np.arange(6)])
    signal = cnv_signal(m)
    assert signal[0] > signal[1] > signal[2]
    corr = cnv_correlation(m, ["t1", "t1", "t1"], signal, top_fraction=0.34)
    assert corr[1] == pytest.approx(1.0)  # same shape as the strongest cell
    assert corr[2] < 0.5


def test_call_malignant_uses_held_out_quantile():
    signal = np.array([0.1, 0.2, 0.3, 5.0])
    corr = np.array([0.0, 0.1, 0.2, 0.9])
    held_out = np.array([True, True, True, False])
    calls, s_thr, c_thr = call_malignant(signal, corr, held_out, quantile=1.0)
    assert calls.tolist() == [False, False, False, True]
    assert s_thr == pytest.approx(0.3) and c_thr == pytest.approx(0.2)


def test_call_malignant_needs_held_out_cells():
    with pytest.raises(ValueError):
        call_malignant(np.ones(2), np.ones(2), np.zeros(2, dtype=bool))

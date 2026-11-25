import numpy as np
import pandas as pd
import pytest

from melscrna.programs import centre_within, derive_program, score_program, split_half


def test_split_half_is_balanced_within_groups():
    groups = pd.Series(["a"] * 10 + ["b"] * 6)
    mask = split_half(groups, seed=0)
    assert mask[groups == "a"].sum() == 5
    assert mask[groups == "b"].sum() == 3


def test_centre_within_removes_group_means():
    data = pd.DataFrame({"g": [1.0, 3.0, 10.0, 14.0]})
    out = centre_within(data, pd.Series(["a", "a", "b", "b"]))
    assert out["g"].tolist() == [-1.0, 1.0, -2.0, 2.0]


def test_derive_program_picks_correlated_genes_and_excludes_anchor():
    rng = np.random.default_rng(0)
    anchor = rng.normal(size=200)
    data = pd.DataFrame(
        {
            "MITF": anchor,
            "TYR": anchor + rng.normal(scale=0.1, size=200),
            "PMEL": anchor + rng.normal(scale=0.5, size=200),
            "NOISE": rng.normal(size=200),
            "AXL": -anchor,
        }
    )
    program = derive_program(data, "MITF", n_genes=2, exclude=["AXL"])
    assert program == ["TYR", "PMEL"]


def test_score_program_is_mean_z_score():
    ref = pd.DataFrame({"a": [0.0, 2.0], "b": [10.0, 20.0]})
    data = pd.DataFrame({"a": [2.0], "b": [10.0]})
    sd_a, sd_b = ref["a"].std(), ref["b"].std()
    expected = ((2.0 - 1.0) / sd_a + (10.0 - 15.0) / sd_b) / 2
    assert score_program(data, ["a", "b"], ref).iloc[0] == pytest.approx(expected)

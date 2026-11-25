"""Phase 6: MITF-high vs AXL-high programmes in CNV-called malignant cells.

Programmes are learned on half the malignant cells (within-tumour correlations) and every
test reported below uses only the other, held-out half, unless stated.
"""

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from melscrna.programs import centre_within, derive_program, score_program, split_half

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD = ROOT / "data" / "processed" / "tirosh_cnv.h5ad"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
NORMAL_TYPES = ["T cell", "B cell", "Plasma cell", "Macrophage", "pDC", "Endothelial", "CAF"]
N_GENES = 50
MIN_CELLS = 10  # tumours with fewer held-out malignant cells are not tested on their own
N_SPLITS = 50  # repeat the learn/test split to see how much results depend on it (D24)


def one_split(expr: pd.DataFrame, tumour: pd.Series, seed: int) -> dict:
    """Learn programmes on one random half, score the other half; return the key numbers."""
    train = split_half(tumour, seed=seed)
    centred = centre_within(expr[train], tumour[train])
    mitf = derive_program(centred, "MITF", N_GENES, exclude=["AXL"])
    axl = derive_program(centred, "AXL", N_GENES, exclude=["MITF"])
    m, a = score_program(expr, mitf, expr), score_program(expr, axl, expr)
    held = ~train
    within = centre_within(pd.DataFrame({"M": m[held], "A": a[held]}), tumour[held])
    frac = ((a - m) > 0).groupby(tumour).mean()
    return {
        "seed": seed,
        "pooled_rho": spearmanr(m[held], a[held]).statistic,
        "within_rho": spearmanr(within["M"], within["A"]).statistic,
        "within_p": spearmanr(within["M"], within["A"]).pvalue,
        **{f"frac_AXL_{t}": v for t, v in frac.items()},
    }


if __name__ == "__main__":
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(IN_H5AD)
    keep = (adata.obs["cnv_call"] == "malignant") & ~adata.obs["cell_type"].isin(NORMAL_TYPES)
    mal = adata[keep.to_numpy()]
    expr = pd.DataFrame(np.asarray(mal.X), index=mal.obs_names, columns=mal.var_names)
    expr = expr.loc[:, (expr > 0).mean() > 0.1]  # genes detected in >10% of malignant cells
    tumour = mal.obs["tumour"].astype(str)
    print(
        f"Malignant cells (CNV call, outside immune/stromal clusters): {len(expr):,} "
        f"from {tumour.nunique()} tumours; genes tested: {expr.shape[1]:,}"
    )

    train = split_half(tumour, seed=0)
    centred_train = centre_within(expr[train], tumour[train])
    mitf_genes = derive_program(centred_train, "MITF", N_GENES, exclude=["AXL"])
    axl_genes = derive_program(centred_train, "AXL", N_GENES, exclude=["MITF"])
    pd.DataFrame({"MITF_programme": mitf_genes, "AXL_programme": axl_genes}).to_csv(
        TABLES / "programme_genes.tsv", sep="\t", index=False
    )
    print(f"\nMITF programme (top 15 of {N_GENES}): {', '.join(mitf_genes[:15])}")
    print(f"AXL programme  (top 15 of {N_GENES}): {', '.join(axl_genes[:15])}")
    print(f"Genes in both programmes: {len(set(mitf_genes) & set(axl_genes))}")

    scores = pd.DataFrame(
        {
            "tumour": tumour,
            "leiden": mal.obs["leiden"].astype(str),
            "held_out": ~train,
            "MITF_score": score_program(expr, mitf_genes, expr),
            "AXL_score": score_program(expr, axl_genes, expr),
        }
    )
    scores["AXL_minus_MITF"] = scores["AXL_score"] - scores["MITF_score"]
    scores.to_csv(TABLES / "programme_scores.tsv", sep="\t")

    test = scores[scores["held_out"]]
    rho_pooled = spearmanr(test["MITF_score"], test["AXL_score"])
    within = centre_within(test[["MITF_score", "AXL_score"]], test["tumour"])
    rho_within = spearmanr(within["MITF_score"], within["AXL_score"])
    print(f"\nHeld-out cells: {len(test)}")
    print(
        f"MITF vs AXL programme, all held-out cells: Spearman rho = "
        f"{rho_pooled.statistic:.2f} (p = {rho_pooled.pvalue:.1e})"
    )
    print(
        f"Within tumours (tumour means removed):     Spearman rho = "
        f"{rho_within.statistic:.2f} (p = {rho_within.pvalue:.1e})"
    )

    rows = []
    for t, grp in scores.groupby("tumour"):
        held = grp[grp["held_out"]]
        row = {
            "tumour": t,
            "n_malignant": len(grp),
            "median_MITF_score": grp["MITF_score"].median(),
            "median_AXL_score": grp["AXL_score"].median(),
            "frac_AXL_dominant": (grp["AXL_minus_MITF"] > 0).mean(),
            "n_held_out": len(held),
        }
        if len(held) >= MIN_CELLS:
            r = spearmanr(held["MITF_score"], held["AXL_score"])
            row.update(within_rho=r.statistic, within_p=r.pvalue)
        rows.append(row)
    per_tumour = pd.DataFrame(rows).set_index("tumour").sort_values("frac_AXL_dominant")
    per_tumour.round(3).to_csv(TABLES / "programmes_per_tumour.tsv", sep="\t")
    with pd.option_context("display.width", 200):
        print("\nPer tumour (all malignant cells; within_rho uses held-out cells only):")
        print(per_tumour.round(2).to_string())

    big = per_tumour[per_tumour["n_malignant"] >= MIN_CELLS]
    r_t = spearmanr(big["median_MITF_score"], big["median_AXL_score"])
    print(
        f"\nBetween tumours ({len(big)} with >= {MIN_CELLS} malignant cells): Spearman rho "
        f"of median MITF vs median AXL score = {r_t.statistic:.2f} (p = {r_t.pvalue:.2f})"
    )

    by_cluster = scores.groupby("leiden")["AXL_minus_MITF"].agg(
        n="count", median="median", frac_AXL_dominant=lambda v: (v > 0).mean()
    )
    print("\nBy Phase 3 cluster:")
    print(by_cluster.round(2).to_string())

    # The numbers above come from one random split (seed 0). Repeat with other splits: a
    # result worth reporting should not depend on which half of the cells was used.
    print(f"\nRepeating with {N_SPLITS} random splits (takes a few minutes)...")
    stab = pd.DataFrame([one_split(expr, tumour, s) for s in range(N_SPLITS)])
    stab.round(4).to_csv(TABLES / "programme_split_stability.tsv", sep="\t", index=False)
    for col in ["pooled_rho", "within_rho"]:
        lo, med, hi = np.percentile(stab[col], [5, 50, 95])
        print(f"{col}: median {med:.2f}, 90% of splits between {lo:.2f} and {hi:.2f}")
    print(
        f"within_rho < 0 in {(stab['within_rho'] < 0).mean():.0%} of splits; "
        f"p < 0.05 in {(stab['within_p'] < 0.05).mean():.0%}"
    )
    frac_cols = [c for c in stab.columns if c.startswith("frac_AXL_")]
    frac = stab[frac_cols].agg(["min", "median", "max"]).T
    frac.index = [c.removeprefix("frac_AXL_") for c in frac.index]
    print("\nFraction AXL-dominant per tumour across splits:")
    print(frac.sort_values("median").round(2).to_string())

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    order = per_tumour.index.tolist()
    cmap = plt.get_cmap("tab20")
    colour = {t: cmap(i % 20) for i, t in enumerate(sorted(order))}
    for i, t in enumerate(order):
        sub = scores[scores["tumour"] == t]
        axes[1].scatter(
            np.full(len(sub), i) + np.random.default_rng(i).uniform(-0.25, 0.25, len(sub)),
            sub["AXL_minus_MITF"],
            s=5,
            alpha=0.5,
        )
    axes[1].axhline(0, color="black", lw=1, ls="--")
    axes[1].set_xticks(range(len(order)), order, rotation=60)
    axes[1].set(ylabel="AXL score minus MITF score", title="Per tumour (above 0 = AXL-dominant)")
    cmap = plt.get_cmap("tab20")
    for i, (t, sub) in enumerate(test.groupby("tumour")):
        axes[0].scatter(
            sub["MITF_score"], sub["AXL_score"], s=6, alpha=0.6, color=cmap(i % 20), label=t
        )
    axes[0].set(
        xlabel="MITF programme score",
        ylabel="AXL programme score",
        title=f"Held-out malignant cells (rho = {rho_pooled.statistic:.2f})",
    )
    axes[0].legend(fontsize=6, ncol=2, frameon=False, markerscale=2)
    fig.tight_layout()
    fig.savefig(FIGURES / "mitf_axl_programmes.png", dpi=150)

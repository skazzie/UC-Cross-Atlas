"""Figure templates for the UC-Cross-Atlas paper, parameterized so they
fill in the moment HB tables arrive.

Each function takes a long-format DataFrame in the schema produced by
code/06_concordance/compute_concordance.py (and equivalents from
08_cross_method, 09_cross_gwas) and returns a matplotlib Figure ready
for saving. Per-atlas ranking plots also accept the per-(atlas, gwas,
tier) scDRS/seismic outputs directly.

Tests for these templates run against the synthetic fixtures already
committed in code/06_concordance/test_compute_concordance.py — so the
templates have green coverage before any HB output arrives.

Functions:
- concordance_heatmap(df_long, value_col, ...) — square cross-atlas
  rho/kappa/jaccard heatmap. The headline 3x3 (or 5x5 for the broad-
  atlas comparator panel) is this function called once on the right
  result table.
- mhc_sensitivity_panel(df_paired_long, ...) — paired plot showing
  per-cell-type score under MHC-excluded vs MHC-included for one
  (atlas, method, tier, gwas).
- method_sensitivity_panel(df_pair_long, ...) — same shape but
  pairing scDRS vs seismic per cell type.
- per_atlas_ranking_plot(df_scores, ...) — single-atlas waterfall
  ranking with FDR bars.
- control_panel(df_traits_long, ...) — side-by-side UC vs SCZ vs
  Height per cell type, matching the rev-7 diagnostic-doc layout.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---- Common style + helpers ----------------------------------------------

def _setup_axes(ax, title: str | None, xlabel: str = "", ylabel: str = "") -> None:
    if title:
        ax.set_title(title, fontsize=10, loc="left")
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=9)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ---- 1. Cross-atlas concordance heatmap -----------------------------------

def concordance_heatmap(
    df_long: pd.DataFrame,
    value_col: str = "spearman_rho",
    atlases: Sequence[str] | None = None,
    title: str | None = None,
    cmap: str = "RdBu_r",
    vmin: float = -1.0, vmax: float = 1.0,
    annotate: bool = True,
):
    """Square N x N cross-atlas heatmap.

    df_long expected columns: atlas_a, atlas_b, <value_col>. Lower
    triangle filled (upper mirror auto). Diagonal = 1 if value_col is a
    correlation metric, else 0.
    """
    import matplotlib.pyplot as plt

    if atlases is None:
        atlases = sorted(set(df_long["atlas_a"]) | set(df_long["atlas_b"]))
    n = len(atlases)
    M = np.full((n, n), np.nan)
    diag = 1.0 if value_col in ("spearman_rho", "jaccard_top5",
                                  "jaccard_top10", "kappa") else 0.0
    for i in range(n):
        M[i, i] = diag
    idx = {a: i for i, a in enumerate(atlases)}
    for _, row in df_long.iterrows():
        a, b = row["atlas_a"], row["atlas_b"]
        if a not in idx or b not in idx:
            continue
        i, j = idx[a], idx[b]
        v = float(row[value_col]) if pd.notna(row[value_col]) else np.nan
        M[i, j] = v
        M[j, i] = v

    fig, ax = plt.subplots(figsize=(0.7 * n + 1.5, 0.6 * n + 1.0))
    im = ax.imshow(M, cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(atlases, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(atlases, fontsize=8)
    _setup_axes(ax, title or f"Cross-atlas {value_col}")
    if annotate:
        for i in range(n):
            for j in range(n):
                if np.isnan(M[i, j]):
                    continue
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center",
                        fontsize=7, color="black" if abs(M[i, j]) < 0.6 else "white")
    fig.colorbar(im, ax=ax, shrink=0.7, label=value_col)
    fig.tight_layout()
    return fig


# ---- 2. MHC sensitivity paired plot ---------------------------------------

def mhc_sensitivity_panel(
    df_paired_long: pd.DataFrame,
    *,
    score_col: str = "score",
    mhc_col: str = "mhc",
    cell_type_col: str = "cell_type",
    atlases: Sequence[str] | None = None,
    title: str | None = None,
):
    """Per-cell-type paired MHC-excl vs MHC-incl plot.

    df_paired_long columns: atlas, cell_type, mhc ('excl' or 'incl'),
    score. Faceted by atlas; cells sorted by MHC-excl score.
    """
    import matplotlib.pyplot as plt

    if atlases is None:
        atlases = sorted(df_paired_long["atlas"].unique())
    nrows = len(atlases)
    fig, axes = plt.subplots(nrows=nrows, figsize=(8, 2.5 * nrows),
                              squeeze=False)
    for k, atlas in enumerate(atlases):
        ax = axes[k, 0]
        sub = df_paired_long[df_paired_long["atlas"] == atlas]
        wide = sub.pivot_table(
            index=cell_type_col, columns=mhc_col, values=score_col,
            aggfunc="first",
        )
        if {"excl", "incl"}.issubset(wide.columns):
            wide = wide.sort_values("excl", ascending=False)
            x = np.arange(len(wide))
            ax.bar(x - 0.2, wide["excl"], width=0.4, label="MHC-excl",
                    color="#7fa3c4")
            ax.bar(x + 0.2, wide["incl"], width=0.4, label="MHC-incl",
                    color="#c47f7f")
            ax.set_xticks(x)
            ax.set_xticklabels(wide.index, rotation=45, ha="right", fontsize=7)
            ax.axhline(0, color="black", linewidth=0.5)
            ax.legend(fontsize=8, loc="best")
        _setup_axes(ax, f"{atlas} — MHC sensitivity",
                    ylabel=score_col)
    if title:
        fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    return fig


# ---- 3. Method sensitivity panel ------------------------------------------

def method_sensitivity_panel(
    df_methods_long: pd.DataFrame,
    *,
    score_col: str = "score",
    method_col: str = "method",
    cell_type_col: str = "cell_type",
    atlases: Sequence[str] | None = None,
    title: str | None = None,
):
    """scDRS vs seismic per cell type, faceted by atlas."""
    import matplotlib.pyplot as plt

    if atlases is None:
        atlases = sorted(df_methods_long["atlas"].unique())
    nrows = len(atlases)
    fig, axes = plt.subplots(nrows=nrows, figsize=(8, 2.5 * nrows),
                              squeeze=False)
    for k, atlas in enumerate(atlases):
        ax = axes[k, 0]
        sub = df_methods_long[df_methods_long["atlas"] == atlas]
        wide = sub.pivot_table(
            index=cell_type_col, columns=method_col, values=score_col,
            aggfunc="first",
        )
        if {"scdrs", "seismic"}.issubset(wide.columns):
            wide = wide.sort_values("scdrs", ascending=False)
            x = np.arange(len(wide))
            ax.bar(x - 0.2, wide["scdrs"], width=0.4, label="scDRS",
                    color="#557ab0")
            ax.bar(x + 0.2, wide["seismic"], width=0.4, label="seismicGWAS",
                    color="#b07a55")
            ax.set_xticks(x)
            ax.set_xticklabels(wide.index, rotation=45, ha="right", fontsize=7)
            ax.axhline(0, color="black", linewidth=0.5)
            ax.legend(fontsize=8, loc="best")
        _setup_axes(ax, f"{atlas} — method sensitivity",
                    ylabel=score_col)
    if title:
        fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    return fig


# ---- 4. Per-atlas ranking waterfall --------------------------------------

def per_atlas_ranking_plot(
    df_scores: pd.DataFrame,
    *,
    score_col: str = "score",
    fdr_col: str | None = "fdr",
    cell_type_col: str = "cell_type",
    title: str | None = None,
    fdr_threshold: float = 0.05,
):
    """Single-atlas ranking waterfall. Highlights cell types passing FDR."""
    import matplotlib.pyplot as plt

    df = df_scores.copy()
    df = df.sort_values(score_col, ascending=True)
    x = np.arange(len(df))
    colors = ["#7c7c7c"] * len(df)
    if fdr_col and fdr_col in df.columns:
        for i, q in enumerate(df[fdr_col]):
            if pd.notna(q) and q < fdr_threshold:
                colors[i] = "#c33"
    fig, ax = plt.subplots(figsize=(0.4 * len(df) + 1.5, 4))
    ax.barh(x, df[score_col], color=colors)
    ax.set_yticks(x)
    ax.set_yticklabels(df[cell_type_col], fontsize=8)
    ax.axvline(0, color="black", linewidth=0.5)
    _setup_axes(ax, title or "Per-atlas ranking", xlabel=score_col)
    fig.tight_layout()
    return fig


# ---- 5. Control side-by-side panel (UC vs SCZ vs Height) ------------------

def control_panel(
    df_traits_long: pd.DataFrame,
    *,
    score_col: str = "score",
    trait_col: str = "trait",
    cell_type_col: str = "cell_type",
    traits: Sequence[str] = ("UC", "SCZ", "Height"),
    cmap: str = "RdBu_r",
    vmin: float | None = None, vmax: float | None = None,
    title: str | None = None,
):
    """Side-by-side trait comparison heatmap (rows = cell types,
    columns = traits). Matches the rev-7 diagnostic-doc layout."""
    import matplotlib.pyplot as plt

    wide = df_traits_long.pivot_table(
        index=cell_type_col, columns=trait_col, values=score_col,
        aggfunc="first",
    )
    wide = wide[[t for t in traits if t in wide.columns]]
    if "UC" in wide.columns:
        wide = wide.sort_values("UC", ascending=False)
    if vmin is None:
        vmin = -np.nanmax(np.abs(wide.values))
    if vmax is None:
        vmax = np.nanmax(np.abs(wide.values))
    fig, ax = plt.subplots(figsize=(0.6 * len(wide.columns) + 1.2,
                                      0.4 * len(wide) + 1.0))
    im = ax.imshow(wide.values, cmap=cmap, vmin=vmin, vmax=vmax,
                    aspect="auto")
    ax.set_xticks(range(len(wide.columns)))
    ax.set_xticklabels(wide.columns, fontsize=9)
    ax.set_yticks(range(len(wide)))
    ax.set_yticklabels(wide.index, fontsize=8)
    _setup_axes(ax, title or "Control panel")
    for i in range(len(wide)):
        for j in range(len(wide.columns)):
            v = wide.values[i, j]
            if pd.isna(v):
                continue
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                    color="black" if abs(v) < 0.6 * max(abs(vmin), abs(vmax)) else "white")
    fig.colorbar(im, ax=ax, shrink=0.6, label=score_col)
    fig.tight_layout()
    return fig


# ---- CLI: smoke-render every template on Garrido data --------------------

def main() -> int:
    """Render every template against the current Garrido outputs.

    Produces results/figures/*.png so the laptop has a visual confirmation
    each template runs cleanly before HB delivery.
    """
    import matplotlib
    matplotlib.use("Agg")

    out_dir = Path(__file__).resolve().parents[2] / "results" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    repo = Path(__file__).resolve().parents[2]

    # 1. Cross-atlas heatmap, from the prototype CSV.
    p = repo / "results/concordance/cross_atlas_pairwise_prototype.csv"
    if p.exists():
        df = pd.read_csv(p)
        fig = concordance_heatmap(df, value_col="spearman_rho",
                                   title="Cross-atlas Spearman rho (Garrido real + 4 synthetic)")
        out = out_dir / "fig_cross_atlas_heatmap.png"
        fig.savefig(out, dpi=150)
        print(f"  wrote {out}")

    # 2 + 3 + 5 need data we'll cobble from the diagnostic outputs.
    rows_traits = []
    for trait, path in [
        ("UC",  repo / "results/scdrs/garrido_delange_seed42/UC.scdrs_group.cell_type_broad"),
        ("SCZ", repo / "results/scdrs/garrido_scz/SCZ.scdrs_group.cell_type_broad"),
        ("Height", repo / "results/scdrs/garrido_height/HEIGHT.scdrs_group.cell_type_broad"),
    ]:
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            rows_traits.append({
                "cell_type": r["group"], "trait": trait,
                "score": r["assoc_mcz"],
            })
    if rows_traits:
        df_traits = pd.DataFrame(rows_traits)
        fig = control_panel(df_traits,
                             title="Garrido controls panel (scDRS, MHC-excl)")
        out = out_dir / "fig_control_panel.png"
        fig.savefig(out, dpi=150)
        print(f"  wrote {out}")

    # 2 — MHC sensitivity (UC excl vs UC incl).
    rows_mhc = []
    for mhc, path in [
        ("excl", repo / "results/scdrs/garrido_delange_seed42/UC.scdrs_group.cell_type_broad"),
        ("incl", repo / "results/scdrs/garrido_delange_mhc/UC_MHC.scdrs_group.cell_type_broad"),
    ]:
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            rows_mhc.append({
                "atlas": "garrido", "cell_type": r["group"], "mhc": mhc,
                "score": r["assoc_mcz"],
            })
    if rows_mhc:
        df_mhc = pd.DataFrame(rows_mhc)
        fig = mhc_sensitivity_panel(df_mhc,
                                     title="Garrido UC MHC sensitivity (scDRS)")
        out = out_dir / "fig_mhc_sensitivity.png"
        fig.savefig(out, dpi=150)
        print(f"  wrote {out}")

    # 3 — method sensitivity (scdrs vs seismic, UC MHC-excl).
    p_sd = repo / "results/scdrs/garrido_delange_seed42/UC.scdrs_group.cell_type_broad"
    p_ss = repo / "results/seismic/garrido_delange_broad.tsv"
    if p_sd.exists() and p_ss.exists():
        sd = pd.read_csv(p_sd, sep="\t")
        ss = pd.read_csv(p_ss, sep="\t")
        rows_method = []
        for _, r in sd.iterrows():
            rows_method.append({"atlas": "garrido", "cell_type": r["group"],
                                "method": "scdrs", "score": r["assoc_mcz"]})
        for _, r in ss.iterrows():
            rows_method.append({"atlas": "garrido", "cell_type": r["cell_type"],
                                "method": "seismic", "score": r["score"]})
        df_method = pd.DataFrame(rows_method)
        fig = method_sensitivity_panel(df_method,
                                        title="Garrido UC method sensitivity")
        out = out_dir / "fig_method_sensitivity.png"
        fig.savefig(out, dpi=150)
        print(f"  wrote {out}")

    # 4 — per-atlas ranking (Garrido UC).
    if p_sd.exists():
        sd = pd.read_csv(p_sd, sep="\t")
        from statsmodels.stats.multitest import multipletests
        _, fdr, _, _ = multipletests(sd["assoc_mcp"], method="fdr_bh")
        df_rank = pd.DataFrame({
            "cell_type": sd["group"],
            "score": sd["assoc_mcz"],
            "fdr": fdr,
        })
        fig = per_atlas_ranking_plot(df_rank,
                                      title="Garrido UC scDRS ranking (MHC-excl)")
        out = out_dir / "fig_per_atlas_ranking.png"
        fig.savefig(out, dpi=150)
        print(f"  wrote {out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

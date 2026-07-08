"""Cross-atlas concordance prototype: Garrido real + 4 synthetic atlases.

Track B.2 of the laptop validation handoff. Validates the cross-atlas
machinery (06_concordance) end-to-end against Garrido's real scoring
output + synthetic stand-ins for the four HB-bound atlases (Smillie,
TAURUS, HCA Gut, Pan-GI). The moment HB delivers real scoring for the
other four atlases, swap out the synthetic generators and rerun — no
code changes.

Inputs the run uses:
  - data/atlases/celltype_crosswalk.tsv (canonical_broad terms per atlas)
  - results/scdrs/garrido_delange_seed42/UC.scdrs_group.cell_type_broad
      (real Garrido scoring)
  - results/seismic/garrido_delange_broad.tsv
  - Synthetic stand-ins for the other four atlases, generated from
    Garrido's scoring with controlled correlation injected per atlas
    so the cross-atlas concordance has a known truth.

Outputs:
  results/concordance/cross_atlas_pairwise_prototype.csv

The synthetic generator uses a fixed seed so the test is deterministic.
The atlas-pair Spearman rho should fall in the rho_target band per
construction; this proves the full PATH:ATLAS:METHOD:GWAS:TIER plumbing
works on real outputs as well as fixtures.

Usage:
  py.exe code/06_concordance/prototype_cross_atlas.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
sys.path.insert(0, str(_HERE))
from metrics import concordance  # noqa: E402


# Atlases participating in the cross-atlas comparison. Real = Garrido;
# synthetic = the other four. Per-atlas rho_target controls how strongly
# the synthetic scoring tracks Garrido. The locked v1 UC trio is
# Smillie + Garrido + TAURUS (HCA Gut + Pan-GI are broad-atlas
# comparators, scope 10/11).
ATLASES = ("garrido", "smillie", "taurus", "hca_gut", "pangi")
SYNTH_RHO_TARGETS = {
    # locked v1 UC trio — strong agreement expected
    "smillie":  0.75,
    "taurus":   0.65,
    # broad-atlas comparators — moderate agreement
    "hca_gut":  0.50,
    "pangi":    0.55,
}


def _load_garrido_real_broad() -> pd.DataFrame:
    """Load Garrido's real scDRS broad-tier output as the anchor."""
    p = _REPO / "results/scdrs/garrido_delange_seed42/UC.scdrs_group.cell_type_broad"
    if not p.exists():
        raise SystemExit(
            f"Garrido real scDRS output missing: {p}\n"
            "Run the Track [1] scDRS dry run first."
        )
    df = pd.read_csv(p, sep="\t")
    df = df.rename(columns={
        "group": "cell_type", "assoc_mcz": "score",
        "assoc_mcp": "pvalue", "n_cell": "n_cells",
    })
    from statsmodels.stats.multitest import multipletests
    _, fdr, _, _ = multipletests(df["pvalue"], method="fdr_bh")
    df["fdr"] = fdr
    return df[["cell_type", "score", "pvalue", "fdr", "n_cells"]].copy()


def _make_synthetic(atlas: str, anchor: pd.DataFrame,
                    rho_target: float, seed: int) -> pd.DataFrame:
    """Generate synthetic per-atlas scoring tracking the anchor at rho ≈ target.

    All atlases share the canonical_broad cell-type vocabulary, so the
    synthetic atlas's cell types are the same 15 _BROAD_VOCAB terms.
    Scores are linear noise on the anchor's z; n_cells perturbed ±50%
    (some cells drop below the 50-cell threshold per the locked v1
    min-cells filter — that's deliberate, tests the filter path).
    """
    rng = np.random.default_rng(seed)
    anchor_z = anchor["score"].to_numpy()
    noise_scale = np.sqrt(max(1e-9, 1 - rho_target ** 2))
    synth_z = rho_target * anchor_z + rng.normal(scale=noise_scale,
                                                  size=len(anchor_z)) * anchor_z.std()
    # P-values from |z| under a fake N(0,1) null, then FDR.
    from scipy.stats import norm
    pvals = 2 * (1 - norm.cdf(np.abs(synth_z)))
    from statsmodels.stats.multitest import multipletests
    _, fdr, _, _ = multipletests(pvals, method="fdr_bh")
    n_cells = (anchor["n_cells"].to_numpy() *
               rng.uniform(0.5, 1.5, size=len(anchor))).astype(int)
    return pd.DataFrame({
        "cell_type": anchor["cell_type"].to_numpy(),
        "score": synth_z,
        "pvalue": pvals,
        "fdr": fdr,
        "n_cells": n_cells,
    })


def _pairwise_concordance(per_atlas: dict[str, pd.DataFrame]) -> list[dict]:
    """Run code/06_concordance/metrics.concordance() on every atlas pair."""
    rows: list[dict] = []
    atlas_keys = sorted(per_atlas)
    for i, a in enumerate(atlas_keys):
        for b in atlas_keys[i + 1:]:
            df_a, df_b = per_atlas[a], per_atlas[b]
            shared = sorted(set(df_a["cell_type"]) & set(df_b["cell_type"]))
            sub_a = df_a[df_a["cell_type"].isin(shared)].set_index("cell_type")
            sub_b = df_b[df_b["cell_type"].isin(shared)].set_index("cell_type")
            cr = concordance(
                scores_a=sub_a["score"].to_dict(),
                scores_b=sub_b["score"].to_dict(),
                qvals_a=sub_a["fdr"].to_dict(),
                qvals_b=sub_b["fdr"].to_dict(),
                cell_counts_a=sub_a["n_cells"].astype(int).to_dict(),
                cell_counts_b=sub_b["n_cells"].astype(int).to_dict(),
                min_cells=50, n_bootstrap=1000, seed=42,
                is_fine_tier=False, larger_is_stronger=True,
            )
            rows.append({
                "atlas_a": a, "atlas_b": b,
                "spearman_rho": cr.spearman_rho,
                "ci_lo": cr.spearman_ci_lo, "ci_hi": cr.spearman_ci_hi,
                "jaccard_top5": cr.jaccard_top5,
                "jaccard_top10": cr.jaccard_top10,
                "kappa": cr.kappa,
                "kappa_threshold": cr.kappa_threshold,
                "kappa_saturated": cr.kappa_threshold_used_due_to_saturation,
                "n_common": cr.n_common,
                "n_excluded_low_count": cr.excluded_low_count,
                "synth_status": (
                    "real-real" if (a == "garrido" and b == "garrido")
                    else "real-synth" if (a == "garrido" or b == "garrido")
                    else "synth-synth"
                ),
            })
    return rows


def main() -> int:
    print("[prototype_cross_atlas] loading Garrido real scoring", flush=True)
    garrido = _load_garrido_real_broad()
    print(f"  Garrido: {len(garrido)} cell types, anchor z range "
          f"[{garrido['score'].min():.2f}, {garrido['score'].max():.2f}]",
          flush=True)

    per_atlas = {"garrido": garrido}
    for i, atlas in enumerate(["smillie", "taurus", "hca_gut", "pangi"]):
        rho = SYNTH_RHO_TARGETS[atlas]
        synth = _make_synthetic(atlas, garrido, rho_target=rho,
                                seed=42 + i * 7)
        per_atlas[atlas] = synth
        print(f"  {atlas:8s}: synthetic, rho_target={rho:.2f}", flush=True)

    rows = _pairwise_concordance(per_atlas)
    out = _REPO / "results/concordance/cross_atlas_pairwise_prototype.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print()
    print(f"[prototype_cross_atlas] wrote {out} ({len(df)} pairs)")
    print(df[["atlas_a", "atlas_b", "spearman_rho", "ci_lo", "ci_hi",
              "n_common", "kappa", "synth_status"]]
          .to_string(index=False, float_format=lambda x: f"{x:7.3f}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

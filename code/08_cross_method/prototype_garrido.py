"""Cross-method concordance prototype on Garrido dry-run outputs.

Track [4] of the laptop validation handoff. Calls the metrics library
on real scDRS + seismic outputs (not the synthetic fixtures the
test suite uses) for Garrido x {de Lange, SCZ, height} x {MHC-excl,
MHC-incl where available}.

The locked v1 run_cross_method.py expects path conventions:
  results/scdrs/<atlas>_<gwas>/cell_type_<tier>/<atlas>_<gwas>.scdrs_group
  results/seismic/<atlas>_<gwas>_<tier>.tsv

Our actual dry-run outputs use a flatter naming:
  results/scdrs/garrido_<gwas>(_mhc|_nocov)/UC(_MHC)?.scdrs_group.cell_type_<tier>
  results/seismic/garrido_<gwas>(_mhc)?_<tier>.tsv

Rather than rename for the dry run, this adapter calls the metrics
directly on the actual files. Outputs a single CSV with one row per
(trait, mhc, tier) comparison.

Usage:
  py.exe code/08_cross_method/prototype_garrido.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
sys.path.insert(0, str(_REPO / "code" / "06_concordance"))
from metrics import concordance  # noqa: E402


def _load_scdrs(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={
        "group": "cell_type", "assoc_mcz": "z_mean",
        "assoc_mcp": "pvalue", "n_cell": "n_cells",
    })
    df["fdr"] = df["pvalue"]  # one trait — scDRS doesn't auto-FDR-correct across cell types here
    from statsmodels.stats.multitest import multipletests
    _, fdr, _, _ = multipletests(df["pvalue"], method="fdr_bh")
    df["fdr"] = fdr
    return df[["cell_type", "z_mean", "pvalue", "fdr", "n_cells"]].copy()


def _load_seismic(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    # seismic ships pvalue + fdr already; coefficient absent in this API.
    # Use score = -log10(pvalue) as the headline magnitude (signed
    # one-sided by construction — slope positive).
    return df[["cell_type", "score", "pvalue", "fdr", "n_cells"]].copy()


def _compare(scdrs_path: Path, seismic_path: Path, label: str) -> dict:
    sd = _load_scdrs(scdrs_path)
    ss = _load_seismic(seismic_path)
    shared = sorted(set(sd["cell_type"]) & set(ss["cell_type"]))
    sd = sd[sd["cell_type"].isin(shared)].set_index("cell_type")
    ss = ss[ss["cell_type"].isin(shared)].set_index("cell_type")
    scores_sd = sd["z_mean"].to_dict()
    scores_ss = ss["score"].to_dict()
    qvals_sd = sd["fdr"].to_dict()
    qvals_ss = ss["fdr"].to_dict()
    counts_sd = sd["n_cells"].astype(int).to_dict()
    counts_ss = ss["n_cells"].astype(int).to_dict()
    cr = concordance(
        scores_a=scores_sd, scores_b=scores_ss,
        qvals_a=qvals_sd, qvals_b=qvals_ss,
        cell_counts_a=counts_sd, cell_counts_b=counts_ss,
        min_cells=50, n_bootstrap=1000, seed=42,
        is_fine_tier=False, larger_is_stronger=True,
    )
    return {
        "label": label,
        "spearman_rho": cr.spearman_rho,
        "ci_lo": cr.spearman_ci_lo,
        "ci_hi": cr.spearman_ci_hi,
        "jaccard_top5": cr.jaccard_top5,
        "jaccard_top10": cr.jaccard_top10,
        "kappa": cr.kappa,
        "kappa_threshold": cr.kappa_threshold,
        "kappa_saturated": cr.kappa_threshold_used_due_to_saturation,
        "n_sig_scdrs": cr.n_sig_a,
        "n_sig_seismic": cr.n_sig_b,
        "n_common": cr.n_common,
    }


def main() -> int:
    cases = [
        # (scdrs_path, seismic_path, label)
        # UC, locked-primary (MHC-excluded, donor cov)
        (_REPO / "results/scdrs/garrido_delange_seed42/UC.scdrs_group.cell_type_broad",
         _REPO / "results/seismic/garrido_delange_broad.tsv",
         "UC_MHCexcl"),
        # UC, MHC-included sensitivity
        (_REPO / "results/scdrs/garrido_delange_mhc/UC_MHC.scdrs_group.cell_type_broad",
         _REPO / "results/seismic/garrido_delange_mhc_broad.tsv",
         "UC_MHCincl"),
        # SCZ, locked-primary — only IF scDRS controls have landed
        (_REPO / "results/scdrs/garrido_scz/SCZ.scdrs_group.cell_type_broad",
         _REPO / "results/seismic/garrido_scz_broad.tsv",
         "SCZ_MHCexcl"),
        # SCZ MHC-included
        (_REPO / "results/scdrs/garrido_scz_mhc/SCZ_MHC.scdrs_group.cell_type_broad",
         _REPO / "results/seismic/garrido_scz_mhc_broad.tsv",
         "SCZ_MHCincl"),
        # Height
        (_REPO / "results/scdrs/garrido_height/HEIGHT.scdrs_group.cell_type_broad",
         _REPO / "results/seismic/garrido_height_broad.tsv",
         "Height_MHCexcl"),
        (_REPO / "results/scdrs/garrido_height_mhc/HEIGHT_MHC.scdrs_group.cell_type_broad",
         _REPO / "results/seismic/garrido_height_mhc_broad.tsv",
         "Height_MHCincl"),
    ]
    rows = []
    for sd_p, ss_p, label in cases:
        if not sd_p.exists():
            print(f"[prototype] SKIP {label}: scDRS output missing ({sd_p.name})",
                  flush=True)
            continue
        if not ss_p.exists():
            print(f"[prototype] SKIP {label}: seismic output missing ({ss_p.name})",
                  flush=True)
            continue
        print(f"[prototype] {label}: {sd_p.name} vs {ss_p.name}", flush=True)
        rows.append(_compare(sd_p, ss_p, label))

    if not rows:
        print("[prototype] no cases available; bailing")
        return 2

    out = _REPO / "results" / "concordance" / "garrido_cross_method_prototype.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print(f"[prototype] wrote {out} ({len(df)} rows)")
    print()
    print(df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

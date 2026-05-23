#!/usr/bin/env python3
"""HCA Gut vs each UC trio atlas concordance, headline + no-Crohn.

Spec: see ./README.md and DECISIONS.md correction 2026-05-20 (6/7).

Loads pre-existing per-atlas method results (scDRS, seismicGWAS) for both
HCA Gut and each UC trio atlas, then computes concordance on the shared
cell-type intersection with the min-50-cells filter.

Does NOT run scDRS or seismicGWAS itself — those are produced by the
respective SLURM wrappers.

For each (uc_atlas, method, tier, gwas) of 3 x 2 x 2 x 2 = 24 combinations,
writes one row to each output TSV (headline vs no-Crohn).
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

LOGGER = logging.getLogger(Path(__file__).stem)

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "code"))
sys.path.insert(0, str(_REPO / "code" / "06_concordance"))

from _shared.constants import (  # noqa: E402
    BOOTSTRAP_N,
    METHODS,
    MIN_CELLS_PER_TYPE,
    SEED,
    TIERS,
    UC_ATLASES,
    UC_GWAS,
)
from _shared.git import git_sha  # noqa: E402
from _shared.result_loading import (  # noqa: E402
    load_scdrs_group,
    load_seismic_results,
    shared_cell_types,
    to_lookup,
)
from metrics import concordance  # noqa: E402


def _load_one(method, scdrs_base, seismic_base, atlas, gwas, tier):
    if method == "scdrs":
        df = load_scdrs_group(scdrs_base, atlas, gwas, tier)
        df = df.rename(columns={"z_mean": "score"})
        return df
    df = load_seismic_results(seismic_base, atlas, gwas, tier)
    from statsmodels.stats.multitest import multipletests
    _, fdr, _, _ = multipletests(df["pvalue"], method="fdr_bh")
    df = df.rename(columns={"coefficient": "score"})
    df["fdr"] = fdr
    return df


def _row(uc_atlas, method, tier, gwas, cr) -> dict:
    return {
        "uc_atlas": uc_atlas,
        "method": method,
        "tier": tier,
        "gwas": gwas,
        "spearman_rho": cr.spearman_rho,
        "ci_lo": cr.spearman_ci_lo,
        "ci_hi": cr.spearman_ci_hi,
        "jaccard_top5": cr.jaccard_top5,
        "jaccard_top10": cr.jaccard_top10,
        "jaccard_top20": cr.jaccard_top20,
        "kappa": cr.kappa,
        "kappa_threshold": cr.kappa_threshold,
        "kappa_saturated": cr.kappa_threshold_used_due_to_saturation,
        "n_sig_hca": cr.n_sig_a,
        "n_sig_uc": cr.n_sig_b,
        "n_common": cr.n_common,
    }


def _run_one_compare(
    method, scdrs_base, seismic_base, hca_results_base,
    uc_atlas, gwas, tier, args,
):
    """Return a single comparison row, or None on skip."""
    try:
        hca_df = _load_one(
            method, hca_results_base, hca_results_base,
            "hca_gut", gwas, tier,
        )
    except FileNotFoundError as exc:
        LOGGER.warning("skip HCA load %s/%s/%s: %s", method, gwas, tier, exc)
        return None
    try:
        uc_df = _load_one(
            method, scdrs_base, seismic_base, uc_atlas, gwas, tier,
        )
    except FileNotFoundError as exc:
        LOGGER.warning("skip UC load %s/%s/%s/%s: %s",
                       uc_atlas, method, gwas, tier, exc)
        return None

    shared = shared_cell_types(hca_df, uc_df, min_cells=args.min_cells)
    if len(shared) < 3:
        LOGGER.warning(
            "skip %s/%s/%s/%s: only %d shared cell types",
            uc_atlas, method, gwas, tier, len(shared),
        )
        return None

    hca_df = hca_df[hca_df["cell_type"].isin(shared)]
    uc_df = uc_df[uc_df["cell_type"].isin(shared)]
    cr = concordance(
        scores_a=to_lookup(hca_df, "cell_type", "score"),
        scores_b=to_lookup(uc_df, "cell_type", "score"),
        qvals_a=to_lookup(hca_df, "cell_type", "fdr"),
        qvals_b=to_lookup(uc_df, "cell_type", "fdr"),
        cell_counts_a=dict(zip(
            hca_df["cell_type"].astype(str),
            hca_df["n_cells"].astype(int),
        )),
        cell_counts_b=dict(zip(
            uc_df["cell_type"].astype(str),
            uc_df["n_cells"].astype(int),
        )),
        min_cells=args.min_cells,
        n_bootstrap=args.bootstrap_n,
        seed=args.seed,
        is_fine_tier=(tier == "fine"),
    )
    return _row(uc_atlas, method, tier, gwas, cr)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--hca-results-base", type=Path, required=True,
                        help="Directory containing HCA Gut scDRS+seismic results "
                             "(headline). Pattern follows the standard layout.")
    parser.add_argument("--hca-no-crohn-base", type=Path,
                        help="Directory with HCA Gut results from the "
                             "no-Crohn-disease sensitivity run.")
    parser.add_argument("--uc-atlases", nargs="+", default=list(UC_ATLASES))
    parser.add_argument("--gwas", nargs="+", default=list(UC_GWAS))
    parser.add_argument("--methods", nargs="+", default=list(METHODS),
                        choices=list(METHODS))
    parser.add_argument("--tiers", nargs="+", default=list(TIERS),
                        choices=["broad", "fine"])
    parser.add_argument("--scdrs-dir", type=Path, required=True)
    parser.add_argument("--seismic-dir", type=Path, required=True)
    parser.add_argument("--out-headline", type=Path, required=True)
    parser.add_argument("--out-no-crohn", type=Path)
    parser.add_argument("--min-cells", type=int, default=MIN_CELLS_PER_TYPE)
    parser.add_argument("--bootstrap-n", type=int, default=BOOTSTRAP_N)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
    )

    import pandas as pd

    t0 = time.time()

    headline_rows = []
    nocrohn_rows = []

    for uc_atlas in args.uc_atlases:
        for method in args.methods:
            for tier in args.tiers:
                for gwas in args.gwas:
                    h_row = _run_one_compare(
                        method, args.scdrs_dir, args.seismic_dir,
                        args.hca_results_base, uc_atlas, gwas, tier, args,
                    )
                    if h_row is not None:
                        headline_rows.append(h_row)
                    if args.out_no_crohn and args.hca_no_crohn_base:
                        nc_row = _run_one_compare(
                            method, args.scdrs_dir, args.seismic_dir,
                            args.hca_no_crohn_base, uc_atlas, gwas, tier, args,
                        )
                        if nc_row is not None:
                            nocrohn_rows.append(nc_row)

    sha = git_sha()
    if not headline_rows:
        LOGGER.error("No headline rows produced.")
        raise SystemExit(2)
    out_h = pd.DataFrame(headline_rows)
    out_h["tool_version"] = "uc-cross-atlas-v1"
    out_h["git_sha"] = sha
    args.out_headline.parent.mkdir(parents=True, exist_ok=True)
    out_h.to_csv(args.out_headline, sep="\t", index=False)
    LOGGER.info("Wrote %s (%d rows)", args.out_headline, len(out_h))

    if args.out_no_crohn and nocrohn_rows:
        out_nc = pd.DataFrame(nocrohn_rows)
        out_nc["tool_version"] = "uc-cross-atlas-v1"
        out_nc["git_sha"] = sha
        args.out_no_crohn.parent.mkdir(parents=True, exist_ok=True)
        out_nc.to_csv(args.out_no_crohn, sep="\t", index=False)
        LOGGER.info("Wrote %s (%d rows)", args.out_no_crohn, len(out_nc))
    elif args.out_no_crohn:
        LOGGER.warning(
            "--out-no-crohn requested but no rows produced; "
            "is --hca-no-crohn-base populated?"
        )

    LOGGER.info("Done in %.1fs", time.time() - t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())

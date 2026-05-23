#!/usr/bin/env python3
"""scDRS vs seismicGWAS within-atlas concordance for (atlas, gwas, tier).

Spec: see ./README.md and code/06_concordance/metrics.py for the locked
metric definitions.

For each (atlas, gwas, tier) of 3 x 2 x 2 = 12 combinations:
  - load scDRS group-analysis output (z_mean, pvalue, fdr, n_cells)
  - load seismicGWAS output (coefficient, se, pvalue, n_cells)
  - filter to shared cell-type intersection with >= 50 cells (from constants)
  - compute Spearman rho + 95% bootstrap CI, top-k Jaccard, Cohen's kappa
    with saturation contingency.

Skips and warns (does not fail) when either input is missing.

Output schema:
  atlas, gwas, tier, spearman_rho, ci_lo, ci_hi,
  jaccard_top5, jaccard_top10, jaccard_top20,
  kappa, kappa_threshold, kappa_saturated,
  n_sig_scdrs, n_sig_seismic, n_common,
  tool_version, git_sha
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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--atlases", nargs="+", default=list(UC_ATLASES))
    parser.add_argument("--gwas", nargs="+", default=list(UC_GWAS))
    parser.add_argument("--tiers", nargs="+", default=list(TIERS),
                        choices=["broad", "fine"])
    parser.add_argument("--scdrs-dir", type=Path, required=True)
    parser.add_argument("--seismic-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
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
    rows = []
    for atlas in args.atlases:
        for gwas in args.gwas:
            for tier in args.tiers:
                try:
                    scdrs_df = load_scdrs_group(
                        args.scdrs_dir, atlas, gwas, tier
                    )
                except FileNotFoundError as exc:
                    LOGGER.warning("skip %s/%s/%s: %s",
                                   atlas, gwas, tier, exc)
                    continue
                try:
                    seis_df = load_seismic_results(
                        args.seismic_dir, atlas, gwas, tier
                    )
                except FileNotFoundError as exc:
                    LOGGER.warning("skip %s/%s/%s: %s",
                                   atlas, gwas, tier, exc)
                    continue

                shared = shared_cell_types(
                    scdrs_df, seis_df, min_cells=args.min_cells
                )
                if len(shared) < 3:
                    LOGGER.warning(
                        "skip %s/%s/%s: only %d shared cell types post-filter",
                        atlas, gwas, tier, len(shared),
                    )
                    continue

                scdrs_df = scdrs_df[scdrs_df["cell_type"].isin(shared)]
                seis_df = seis_df[seis_df["cell_type"].isin(shared)]

                scores_a = to_lookup(scdrs_df, "cell_type", "z_mean")
                scores_b = to_lookup(seis_df, "cell_type", "coefficient")
                qvals_a = to_lookup(scdrs_df, "cell_type", "fdr")
                # seismic outputs raw pvalue; BH-FDR within tool.
                from statsmodels.stats.multitest import multipletests
                _, fdr_b, _, _ = multipletests(
                    seis_df["pvalue"], method="fdr_bh"
                )
                qvals_b = dict(zip(
                    seis_df["cell_type"].astype(str), fdr_b.astype(float)
                ))
                counts_a = dict(zip(
                    scdrs_df["cell_type"].astype(str),
                    scdrs_df["n_cells"].astype(int),
                ))
                counts_b = dict(zip(
                    seis_df["cell_type"].astype(str),
                    seis_df["n_cells"].astype(int),
                ))

                cr = concordance(
                    scores_a=scores_a,
                    scores_b=scores_b,
                    qvals_a=qvals_a,
                    qvals_b=qvals_b,
                    cell_counts_a=counts_a,
                    cell_counts_b=counts_b,
                    min_cells=args.min_cells,
                    n_bootstrap=args.bootstrap_n,
                    seed=args.seed,
                    is_fine_tier=(tier == "fine"),
                )

                rows.append({
                    "atlas": atlas,
                    "gwas": gwas,
                    "tier": tier,
                    "spearman_rho": cr.spearman_rho,
                    "ci_lo": cr.spearman_ci_lo,
                    "ci_hi": cr.spearman_ci_hi,
                    "jaccard_top5": cr.jaccard_top5,
                    "jaccard_top10": cr.jaccard_top10,
                    "jaccard_top20": cr.jaccard_top20,
                    "kappa": cr.kappa,
                    "kappa_threshold": cr.kappa_threshold,
                    "kappa_saturated": cr.kappa_threshold_used_due_to_saturation,
                    "n_sig_scdrs": cr.n_sig_a,
                    "n_sig_seismic": cr.n_sig_b,
                    "n_common": cr.n_common,
                })
                LOGGER.info(
                    "%s/%s/%s: rho=%.3f n_common=%d",
                    atlas, gwas, tier, cr.spearman_rho, cr.n_common,
                )

    if not rows:
        LOGGER.error("No (atlas, gwas, tier) combinations produced output.")
        raise SystemExit(2)

    out = pd.DataFrame(rows)
    out["tool_version"] = "uc-cross-atlas-v1"
    out["git_sha"] = git_sha()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    LOGGER.info("Wrote %s (%d rows) in %.1fs",
                args.out, len(out), time.time() - t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""de Lange vs Liu within-(atlas, method, tier) cross-GWAS concordance.

Spec: see ./README.md.

For each (atlas, method, tier) of 3 x 2 x 2 = 12 combinations:
  - load cell-type-level results under de Lange and under Liu
  - compute Spearman rho + 95% bootstrap CI, top-k Jaccard, Cohen's kappa
    with saturation contingency.

Skips and warns (does not fail) when either input is missing.

Output schema:
  atlas, method, tier, spearman_rho, ci_lo, ci_hi,
  jaccard_top5, jaccard_top10, jaccard_top20,
  kappa, kappa_threshold, kappa_saturated,
  n_sig_delange, n_sig_liu, n_common,
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
    METHODS,
    MIN_CELLS_PER_TYPE,
    SEED,
    TIERS,
    UC_ATLASES,
)
from _shared.git import git_sha  # noqa: E402
from _shared.result_loading import (  # noqa: E402
    load_scdrs_group,
    load_seismic_results,
    shared_cell_types,
    to_lookup,
)
from metrics import concordance  # noqa: E402


def _load_one(method, base_dir, atlas, gwas, tier):
    if method == "scdrs":
        df = load_scdrs_group(base_dir, atlas, gwas, tier)
        return df.rename(columns={"z_mean": "score"}), df["fdr"]
    df = load_seismic_results(base_dir, atlas, gwas, tier)
    from statsmodels.stats.multitest import multipletests
    _, fdr, _, _ = multipletests(df["pvalue"], method="fdr_bh")
    df = df.rename(columns={"coefficient": "score"})
    df["fdr"] = fdr
    return df, df["fdr"]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--atlases", nargs="+", default=list(UC_ATLASES))
    parser.add_argument("--methods", nargs="+", default=list(METHODS),
                        choices=list(METHODS))
    parser.add_argument("--tiers", nargs="+", default=list(TIERS),
                        choices=["broad", "fine"])
    parser.add_argument("--scdrs-dir", type=Path, required=True)
    parser.add_argument("--seismic-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--gwas-a", default="delange")
    parser.add_argument("--gwas-b", default="liu")
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
        for method in args.methods:
            base_dir = args.scdrs_dir if method == "scdrs" else args.seismic_dir
            for tier in args.tiers:
                try:
                    df_a, _ = _load_one(method, base_dir, atlas, args.gwas_a, tier)
                    df_b, _ = _load_one(method, base_dir, atlas, args.gwas_b, tier)
                except FileNotFoundError as exc:
                    LOGGER.warning("skip %s/%s/%s: %s",
                                   atlas, method, tier, exc)
                    continue

                shared = shared_cell_types(df_a, df_b, min_cells=args.min_cells)
                if len(shared) < 3:
                    LOGGER.warning(
                        "skip %s/%s/%s: %d shared cell types",
                        atlas, method, tier, len(shared),
                    )
                    continue
                df_a = df_a[df_a["cell_type"].isin(shared)]
                df_b = df_b[df_b["cell_type"].isin(shared)]

                cr = concordance(
                    scores_a=to_lookup(df_a, "cell_type", "score"),
                    scores_b=to_lookup(df_b, "cell_type", "score"),
                    qvals_a=to_lookup(df_a, "cell_type", "fdr"),
                    qvals_b=to_lookup(df_b, "cell_type", "fdr"),
                    cell_counts_a=dict(zip(
                        df_a["cell_type"].astype(str),
                        df_a["n_cells"].astype(int),
                    )),
                    cell_counts_b=dict(zip(
                        df_b["cell_type"].astype(str),
                        df_b["n_cells"].astype(int),
                    )),
                    min_cells=args.min_cells,
                    n_bootstrap=args.bootstrap_n,
                    seed=args.seed,
                    is_fine_tier=(tier == "fine"),
                )
                rows.append({
                    "atlas": atlas,
                    "method": method,
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
                    "n_sig_delange": cr.n_sig_a,
                    "n_sig_liu": cr.n_sig_b,
                    "n_common": cr.n_common,
                })
                LOGGER.info(
                    "%s/%s/%s: rho=%.3f n_common=%d",
                    atlas, method, tier, cr.spearman_rho, cr.n_common,
                )

    if not rows:
        LOGGER.error("No (atlas, method, tier) combinations produced output.")
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

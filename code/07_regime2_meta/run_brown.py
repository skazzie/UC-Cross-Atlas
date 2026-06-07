#!/usr/bin/env python3
"""Brown's-method meta-analysis across UC-trio atlases (regime 2).

Spec: see ./README.md and DECISIONS.md.

For one (method, tier, gwas), combine per-cell-type p-values across
{smillie, garrido_trigo, taurus} using Brown's method (Kost-McDermott
analytical form), with the pairwise covariance estimated empirically
from null draws.

Implementation note: the spec referenced an `EmpiricalBrownsMethod`
Python package, but no such package exists on PyPI. Implementing the
algorithm directly here is ~30 lines and avoids an unsatisfiable
dependency. The math matches the R `EmpiricalBrownsMethod::kostsMethod()`
implementation (Poole et al. 2016, Bioinformatics) using Brown's (1975)
polynomial approximation for cov(-2 ln p_i, -2 ln p_j) in terms of the
underlying test-statistic correlation.

Output TSV columns:
  cell_type, combined_pval, n_atlases_combined, correlation_fallback, method,
  tier, gwas, tool_version, git_sha
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

from _shared.constants import UC_ATLASES  # noqa: E402
from _shared.git import git_sha  # noqa: E402
from _shared.result_loading import (  # noqa: E402
    load_scdrs_group,
    load_seismic_results,
)


def require_path(p: Path, descr: str) -> None:
    if not p.exists():
        LOGGER.error("Missing input: %s (%s)", p, descr)
        raise SystemExit(2)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--method", required=True, choices=["scdrs", "seismic"])
    parser.add_argument("--tier", required=True, choices=["broad", "fine"])
    parser.add_argument("--gwas", default="delange")
    parser.add_argument("--regime1-dir", type=Path, required=True,
                        help="results/scdrs or results/seismic")
    parser.add_argument("--null-draws-dir", type=Path, required=True,
                        help="results/null_draws")
    parser.add_argument("--atlases", nargs="+", default=list(UC_ATLASES))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
    )

    import numpy as np
    import pandas as pd
    from scipy.stats import chi2

    def _brown_cov_from_corr(r: float) -> float:
        """Brown's (1975) polynomial approximation: cov(-2 ln p_i, -2 ln p_j)
        as a function of the Pearson correlation of the underlying test
        statistics. Valid for normally-distributed test statistics."""
        return 3.263 * r + 0.710 * r * r + 0.027 * r ** 3

    def brown_combine(pvals: np.ndarray, corr: np.ndarray) -> float:
        """Kost-McDermott combined p-value.

        pvals: 1-D array of k input p-values.
        corr:  k x k correlation matrix of the underlying test statistics
               (Pearson). Diagonal == 1.0.
        """
        pvals = np.clip(pvals, 1e-300, 1.0)
        k = len(pvals)
        if k == 1:
            return float(pvals[0])
        T = -2.0 * np.sum(np.log(pvals))
        E = 2.0 * k  # expected value under independence
        var = 4.0 * k  # variance under independence
        for i in range(k):
            for j in range(i + 1, k):
                var += 2.0 * _brown_cov_from_corr(float(corr[i, j]))
        c = var / (2.0 * E)
        f = 2.0 * E * E / var
        return float(chi2.sf(T / c, df=f))

    t0 = time.time()

    # Per-atlas cell-type-level p-values + null tensors.
    per_atlas_p: dict[str, pd.DataFrame] = {}
    per_atlas_nulls: dict[str, dict[str, np.ndarray]] = {}

    for atlas in args.atlases:
        if args.method == "scdrs":
            try:
                df = load_scdrs_group(
                    args.regime1_dir, atlas, args.gwas, args.tier
                )
                df = df.rename(columns={"pvalue": "p"})
            except FileNotFoundError:
                LOGGER.warning(
                    "Missing scDRS group output for %s; skipping atlas.", atlas
                )
                continue
        else:
            try:
                df = load_seismic_results(
                    args.regime1_dir, atlas, args.gwas, args.tier
                )
                df = df.rename(columns={"pvalue": "p"})
            except FileNotFoundError:
                LOGGER.warning(
                    "Missing seismic output for %s; skipping atlas.", atlas
                )
                continue
        per_atlas_p[atlas] = df

        null_path = args.null_draws_dir / f"{atlas}_{args.tier}_nulls.npz"
        if not null_path.exists():
            LOGGER.warning(
                "Missing null tensor for %s at %s; atlas excluded from "
                "correlation estimation.", atlas, null_path,
            )
            continue
        npz = np.load(null_path, allow_pickle=True)
        cts = np.asarray(npz["cell_types"], dtype=object)
        nulls = np.asarray(npz["nulls"], dtype=np.float32)
        per_atlas_nulls[atlas] = {
            str(ct): nulls[i] for i, ct in enumerate(cts)
        }
        LOGGER.info("Loaded %s nulls: %d cell types x %d draws",
                    atlas, len(cts), nulls.shape[1])

    if not per_atlas_p:
        LOGGER.error("No per-atlas results loaded. Aborting.")
        raise SystemExit(2)

    # All cell types appearing in at least one atlas.
    all_cts: set[str] = set()
    for df in per_atlas_p.values():
        all_cts.update(df["cell_type"].astype(str))
    LOGGER.info("Cell types across atlases: %d", len(all_cts))

    # Edge-case fallback: cell types whose null SD is below the 5th
    # percentile in any atlas use the median cross-atlas correlation.
    null_sd_per_atlas: dict[str, dict[str, float]] = {}
    median_corr_per_atlas: dict[str, float] = {}
    for atlas, nulls in per_atlas_nulls.items():
        sds = {ct: float(np.std(v, ddof=1)) for ct, v in nulls.items()}
        null_sd_per_atlas[atlas] = sds

    # Pairwise median correlation across atlas pairs, over all common
    # cell types — used as the substitute when an atlas's null SD is
    # pathologically low for a specific cell type.
    atlas_keys = list(per_atlas_nulls)
    pair_corrs: list[float] = []
    for i, a in enumerate(atlas_keys):
        for b in atlas_keys[i + 1:]:
            shared = set(per_atlas_nulls[a]) & set(per_atlas_nulls[b])
            for ct in shared:
                v_a = per_atlas_nulls[a][ct]
                v_b = per_atlas_nulls[b][ct]
                n = min(len(v_a), len(v_b))
                if n < 10:
                    continue
                r = float(np.corrcoef(v_a[:n], v_b[:n])[0, 1])
                if np.isfinite(r):
                    pair_corrs.append(r)
    median_cross_corr = float(np.median(pair_corrs)) if pair_corrs else 0.0
    LOGGER.info("Median cross-atlas null correlation: %.4f (n=%d pairs)",
                median_cross_corr, len(pair_corrs))

    # 5th-percentile SD per atlas.
    sd_5p: dict[str, float] = {}
    for atlas, sds in null_sd_per_atlas.items():
        vals = np.array(list(sds.values()), dtype=float)
        sd_5p[atlas] = float(np.percentile(vals, 5)) if vals.size else 0.0

    def build_corr_matrix(
        ct: str, atlases_present: list[str]
    ) -> tuple[np.ndarray, bool]:
        """Empirical correlation matrix for the null draws of one cell type.

        Returns (matrix, fallback_used).
        """
        k = len(atlases_present)
        M = np.eye(k, dtype=float)
        fallback = False
        for i in range(k):
            for j in range(i + 1, k):
                a = atlases_present[i]
                b = atlases_present[j]
                v_a = per_atlas_nulls.get(a, {}).get(ct)
                v_b = per_atlas_nulls.get(b, {}).get(ct)
                if v_a is None or v_b is None:
                    M[i, j] = M[j, i] = median_cross_corr
                    fallback = True
                    continue
                sd_a_low = null_sd_per_atlas[a].get(ct, 0) < sd_5p[a]
                sd_b_low = null_sd_per_atlas[b].get(ct, 0) < sd_5p[b]
                if sd_a_low or sd_b_low:
                    M[i, j] = M[j, i] = median_cross_corr
                    fallback = True
                    continue
                n = min(len(v_a), len(v_b))
                r = float(np.corrcoef(v_a[:n], v_b[:n])[0, 1])
                if not np.isfinite(r):
                    r = median_cross_corr
                    fallback = True
                if r <= 0:
                    LOGGER.warning(
                        "%s vs %s for %s: non-positive off-diagonal %.4f",
                        a, b, ct, r,
                    )
                M[i, j] = M[j, i] = r
        return M, fallback

    rows = []
    for ct in sorted(all_cts):
        atlases_with_p = [
            a for a, df in per_atlas_p.items()
            if ct in set(df["cell_type"].astype(str))
        ]
        ps = []
        for a in atlases_with_p:
            df = per_atlas_p[a]
            ps.append(
                float(df.loc[df["cell_type"].astype(str) == ct, "p"].iloc[0])
            )
        n_combined = len(atlases_with_p)
        fallback = False
        if n_combined == 1:
            combined = ps[0]
        elif n_combined >= 2:
            atlases_for_corr = [
                a for a in atlases_with_p if a in per_atlas_nulls
                and ct in per_atlas_nulls[a]
            ]
            if len(atlases_for_corr) < n_combined:
                fallback = True
                atlases_for_corr = atlases_with_p
            corr_matrix, fb = build_corr_matrix(ct, atlases_for_corr)
            fallback = fallback or fb
            combined = brown_combine(np.array(ps, dtype=float), corr_matrix)
        else:
            combined = float("nan")
        rows.append({
            "cell_type": ct,
            "combined_pval": combined,
            "n_atlases_combined": n_combined,
            "correlation_fallback": bool(fallback),
            "atlases": ",".join(atlases_with_p),
        })

    result = pd.DataFrame(rows)
    result["method"] = args.method
    result["tier"] = args.tier
    result["gwas"] = args.gwas
    result["tool_version"] = "brown-kost-mcdermott-inline"
    result["git_sha"] = git_sha()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, sep="\t", index=False)
    LOGGER.info(
        "Wrote %s (%d rows, %d with fallback) in %.1fs",
        args.out, len(result),
        int(result["correlation_fallback"].sum()),
        time.time() - t0,
    )

    # TODO (stretch #1): if correlation matrices look pathological, activate
    # the 300-permutation scDRS fallback per README. Decision deferred to M5.
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Run LDSC --h2 and extract intercept / ratio into a one-row TSV.

Wrapper around `ldsc.py --h2`. Produces:
  - <out>.log         LDSC's full output (written by LDSC itself)
  - <out>.intercept.tsv  one-row TSV with the F10 gate fields
                         (gwas, h2, h2_se, intercept, intercept_se,
                          ratio, ratio_se, lambda_gc, mean_chisq,
                          status, in_band)

`status` is one of {ok, partial, missing, missing-file} from
`parse_log.parse()`. `in_band` is True iff F10 acceptance passes
(intercept <= 1.10 OR ratio <= 0.20).

Runs in the `uc-cross-atlas` Python 3 env; LDSC itself is invoked
via subprocess under `--ldsc-python` (Py 2.7 in
`uc-cross-atlas-ldsc`).
"""

from __future__ import annotations

import argparse
import logging
import shlex
import subprocess
import sys
from pathlib import Path

LOGGER = logging.getLogger(Path(__file__).stem)

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from parse_log import in_band_for_f10, parse  # noqa: E402


def build_ldsc_command(
    ldsc_dir: Path, sumstats: Path, ld_dir: Path, out_prefix: Path,
    *, python_bin: str = "python"
) -> list[str]:
    """Build the argv for LDSC's ldsc.py --h2."""
    return [
        python_bin,
        str(ldsc_dir / "ldsc.py"),
        "--h2", str(sumstats),
        "--ref-ld-chr", str(ld_dir) + "/",
        "--w-ld-chr", str(ld_dir) + "/",
        "--out", str(out_prefix),
    ]


def write_intercept_tsv(out_prefix: Path, gwas: str, parsed: dict) -> Path:
    """Write the one-row intercept TSV that feeds F10 sign-off.

    Always called (even on parse failure) so the result is auditable.
    """
    import pandas as pd

    row = {
        "gwas": gwas,
        "h2": parsed.get("h2"),
        "h2_se": parsed.get("h2_se"),
        "intercept": parsed.get("intercept"),
        "intercept_se": parsed.get("intercept_se"),
        "ratio": parsed.get("ratio"),
        "ratio_se": parsed.get("ratio_se"),
        "lambda_gc": parsed.get("lambda_gc"),
        "mean_chisq": parsed.get("mean_chisq"),
        "status": parsed.get("status"),
        "in_band": in_band_for_f10(
            parsed.get("intercept"), parsed.get("ratio")
        ),
    }
    tsv_path = Path(str(out_prefix) + ".intercept.tsv")
    pd.DataFrame([row]).to_csv(tsv_path, sep="\t", index=False)
    return tsv_path


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--sumstats", type=Path, required=True,
                   help=".sumstats.gz from munge_for_ldsc.py")
    p.add_argument("--ld-dir", type=Path, required=True,
                   help="Directory containing eur_w_ld_chr/ LD scores")
    p.add_argument("--out", type=Path, required=True,
                   help="LDSC --out prefix; writes <out>.log + "
                        "<out>.intercept.tsv")
    p.add_argument("--ldsc-dir", type=Path, required=True,
                   help="Path to cloned bulik/ldsc")
    p.add_argument("--gwas", required=True,
                   help="Short tag for the GWAS (e.g. uc_delange)")
    p.add_argument("--ldsc-python", default="python",
                   help="Python interpreter for LDSC (Py 2.7)")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cmd = build_ldsc_command(
        args.ldsc_dir, args.sumstats, args.ld_dir, args.out,
        python_bin=args.ldsc_python,
    )
    LOGGER.info("LDSC command: %s", " ".join(shlex.quote(c) for c in cmd))

    if args.dry_run:
        LOGGER.info("Dry run; not invoking LDSC.")
        return 0

    rc = subprocess.call(cmd)
    if rc != 0:
        LOGGER.warning("LDSC exited %d; will still attempt to parse .log", rc)

    log_path = Path(str(args.out) + ".log")
    parsed = parse(log_path)
    tsv_path = write_intercept_tsv(args.out, args.gwas, parsed)
    LOGGER.info(
        "Wrote %s; intercept=%s ratio=%s status=%s in_band=%s",
        tsv_path, parsed.get("intercept"), parsed.get("ratio"),
        parsed.get("status"), in_band_for_f10(
            parsed.get("intercept"), parsed.get("ratio")
        ),
    )
    return rc


if __name__ == "__main__":
    sys.exit(main())

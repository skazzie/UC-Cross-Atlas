"""Munge prepare_gwas.py intermediates → LDSC sumstats format.

`code/01_magma/prepare_gwas.py` writes:
  - <prefix>.snp.loc  (SNP, CHR, BP)   no header
  - <prefix>.pval     (SNP, P, N)      with header

LDSC's munge_sumstats.py wants one TSV with SNP, A1, A2, N, P (and
optionally Z or beta + se). We join the two prepare_gwas outputs,
optionally pick up A1/A2 from the original harmonized sumstats if
available, write a merged TSV, then shell out to LDSC's
munge_sumstats.py.

The wrapper is a Python 3 script that runs in the `uc-cross-atlas`
env; LDSC's munge_sumstats.py runs in `uc-cross-atlas-ldsc` (Py2.7).
We invoke LDSC via the `--ldsc-dir` flag pointing at the cloned
bulik/ldsc repo.

A1/A2 handling:
  - If the original sumstats file path is passed via --original-sumstats,
    A1/A2 are pulled from it via a join on SNP id. This is the path
    LDSC prefers; the alleles are required for proper strand
    matching to the LD scores.
  - Otherwise, A1/A2 are written as N/N placeholders and LDSC is
    invoked with --no-alleles. This costs a small amount of SNP
    filtering precision but the heritability/intercept estimates
    survive.
"""

from __future__ import annotations

import argparse
import logging
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

LOGGER = logging.getLogger(Path(__file__).stem)


def build_merged_sumstats(snp_loc: Path, pval: Path,
                          original_sumstats: Path | None,
                          original_a1_col: str | None,
                          original_a2_col: str | None,
                          original_snp_col: str | None) -> "pd.DataFrame":
    """Join .snp.loc + .pval + (optional) original A1/A2 into one DF.

    Returns a DataFrame with columns SNP, CHR, BP, N, P (+ A1, A2 if
    available). Imported here (not at module top) to keep the import
    surface for tests light.
    """
    import pandas as pd

    snp_df = pd.read_csv(snp_loc, sep="\t", header=None,
                         names=["SNP", "CHR", "BP"])
    pval_df = pd.read_csv(pval, sep="\t")
    if "SNP" not in pval_df.columns or "P" not in pval_df.columns:
        raise SystemExit(f"{pval}: expected columns SNP, P, N")
    merged = snp_df.merge(pval_df, on="SNP", how="inner")
    if "N" not in merged.columns:
        raise SystemExit(f"{pval}: no N column; prepare_gwas should have written it")

    if original_sumstats is not None:
        if not (original_a1_col and original_a2_col and original_snp_col):
            raise SystemExit(
                "--original-sumstats requires --original-snp-col, "
                "--original-a1-col, --original-a2-col"
            )
        sep = "\t"
        if original_sumstats.suffix == ".gz" or original_sumstats.name.endswith(".tsv.gz"):
            comp = "gzip"
        else:
            comp = None
        orig = pd.read_csv(
            original_sumstats, sep=sep, compression=comp,
            usecols=[original_snp_col, original_a1_col, original_a2_col],
            low_memory=False,
        )
        orig = orig.rename(columns={
            original_snp_col: "SNP",
            original_a1_col: "A1",
            original_a2_col: "A2",
        })
        merged = merged.merge(orig, on="SNP", how="inner")
        LOGGER.info("After A1/A2 join: %d SNPs", len(merged))
    return merged


def build_ldsc_command(
    ldsc_dir: Path, snplist: Path, in_tsv: Path, out_prefix: Path,
    *, has_alleles: bool, signed_sumstats: str = "P,0"
) -> list[str]:
    """Build the argv for LDSC's munge_sumstats.py.

    `signed_sumstats` defaults to "P,0" (no signed stat — purely from
    p-values; LDSC's intercept estimation tolerates unsigned). For
    studies with beta/OR available, prefer the signed form ("BETA,0").
    """
    cmd = [
        sys.executable,  # placeholder; user overrides with python2.7 in the LDSC env
        str(ldsc_dir / "munge_sumstats.py"),
        "--sumstats", str(in_tsv),
        "--out", str(out_prefix),
        "--merge-alleles", str(snplist),
        "--signed-sumstats", signed_sumstats,
        "--snp", "SNP",
        "--p", "P",
        "--N-col", "N",
    ]
    if has_alleles:
        cmd += ["--a1", "A1", "--a2", "A2"]
    else:
        cmd += ["--no-alleles"]
    return cmd


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--in-snp-loc", type=Path, required=True)
    p.add_argument("--in-pval", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True,
                   help="LDSC --out prefix; LDSC appends .sumstats.gz")
    p.add_argument("--ldsc-dir", type=Path, required=True,
                   help="Path to cloned bulik/ldsc")
    p.add_argument("--snplist", type=Path, required=True,
                   help="Path to w_hm3.snplist")
    p.add_argument("--trait", required=True,
                   help="Short trait tag (UC / SCZ / height); used in logs")
    p.add_argument("--original-sumstats", type=Path, default=None,
                   help="Optional: original harmonized sumstats for A1/A2")
    p.add_argument("--original-snp-col", default=None)
    p.add_argument("--original-a1-col", default=None)
    p.add_argument("--original-a2-col", default=None)
    p.add_argument("--ldsc-python", default="python",
                   help="Python interpreter to run LDSC under "
                        "(should be the uc-cross-atlas-ldsc env's python2.7)")
    p.add_argument("--dry-run", action="store_true",
                   help="Build everything but don't invoke LDSC")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
    )

    LOGGER.info("Merging %s + %s", args.in_snp_loc, args.in_pval)
    merged = build_merged_sumstats(
        args.in_snp_loc, args.in_pval, args.original_sumstats,
        args.original_a1_col, args.original_a2_col, args.original_snp_col,
    )
    LOGGER.info("Merged sumstats: %d rows", len(merged))

    with tempfile.NamedTemporaryFile(
        "w", suffix=".tsv", delete=False, encoding="utf-8"
    ) as tmpf:
        merged.to_csv(tmpf, sep="\t", index=False)
        tmp_path = Path(tmpf.name)
    LOGGER.info("Wrote merged TSV to %s", tmp_path)

    has_alleles = "A1" in merged.columns and "A2" in merged.columns
    args.out.parent.mkdir(parents=True, exist_ok=True)
    cmd = build_ldsc_command(
        args.ldsc_dir, args.snplist, tmp_path, args.out,
        has_alleles=has_alleles,
    )
    cmd[0] = args.ldsc_python
    LOGGER.info("LDSC command: %s", " ".join(shlex.quote(c) for c in cmd))

    if args.dry_run:
        LOGGER.info("Dry run; not invoking LDSC.")
        return 0
    rc = subprocess.call(cmd)
    if rc != 0:
        LOGGER.error("LDSC munge_sumstats exited %d", rc)
        return rc
    LOGGER.info("LDSC munge_sumstats done")
    return 0


if __name__ == "__main__":
    sys.exit(main())

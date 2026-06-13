# 01b_ldsc — LDSC intercept (F10 pre-narrative gate)

**Status.** Code-complete scaffold. Execution is gated on the
Python-2.7 LDSC conda env (separate from the `uc-cross-atlas` env,
because LDSC is Python 2) + WSL/Ubuntu for the LDSC binary. The
wrappers here build the LDSC command lines and run them; LDSC itself
is upstream code.

## Why this module exists (OPEN_FLAGS F10)

The DECISIONS 25 framing: λ_GC at large N is dominated by polygenic
signal, not stratification (Bulik-Sullivan 2015 *Nat Genet*). The
**LDSC intercept** is the correct stratification check. For the
primary UC GWAS (de Lange GCST004133), the intercept must be in-band
before the broad-tier 3×3 heatmap can be written up as biology —
without it, the heatmap is a methods result, not a biology result.

F10 does NOT gate heatmap *generation* (de Lange is the same GWAS
across all five atlases, so any stratification is a shared input
that doesn't differentially distort cross-atlas agreement). It DOES
gate the result paragraph in the M4 manuscript draft.

## Pipeline

```
data/gwas/uc_delange.{snp.loc,pval}    ──┐
data/gwas/yengo_height.{snp.loc,pval}  ──┼──> munge_for_ldsc.py ──> *.sumstats.gz
data/gwas/scz_trubetskoy.{snp.loc,pval} ─┘                              │
                                                                       v
                              eur_w_ld_chr/, w_hm3.snplist  ──> run_ldsc_h2.py
                                                                       │
                                                                       v
                                              results/ldsc/<gwas>.log
                                              (extract intercept, ratio)
```

- `munge_for_ldsc.py` calls LDSC's `munge_sumstats.py` from the
  `.snp.loc` + `.pval` intermediates this repo's `prepare_gwas.py`
  already produces (DECISIONS 24, 26, 28). The wrapper handles the
  N_eff-as-N convention discipline (DECISIONS 23(a), 27(a), 29(a))
  so the LDSC run sees the same effective N MAGMA does.
- `run_ldsc_h2.py` runs `ldsc.py --h2 --ref-ld-chr eur_w_ld_chr/
  --w-ld-chr eur_w_ld_chr/`, scrapes the `.log` for
  `Intercept: X (SE)` and `Ratio: Y (SE)`, writes a per-GWAS results
  TSV ready to plot into OPEN_FLAGS F10 sign-off.

## Acceptance criteria

- de Lange `Intercept` ≤ 1.10 (or `Ratio` ≤ 0.20). Either passes F10
  for the primary GWAS.
- Yengo: report-only (positive control — intercept will be elevated,
  not actionable).
- SCZ: report-only (negative control).
- Liu: gated on the ancestry-LD decision (DECISIONS 27(c)).

## One-time setup (when WSL+Ubuntu is up)

LDSC is Python 2.7 and brittle on modern systems. The cleanest install
is its own conda env under WSL:

```bash
# Inside WSL Ubuntu, from the repo root on /mnt/c/...
cd /mnt/c/Users/muska/UC-Cross-Atlas/code/01b_ldsc
conda env create -f environment_py27.yml
conda activate uc-cross-atlas-ldsc

# Clone bulik/ldsc — pinned to the commit used for v1
git clone https://github.com/bulik/ldsc.git ~/ldsc
cd ~/ldsc && git checkout aa33296abac9569a6422ee6ba7eb4b902422cc74 && cd -

# LD scores — EUR (DECISIONS 14 + this module's README):
# eur_w_ld_chr.tar.bz2, w_hm3.snplist.bz2
# Mirrored at https://alkesgroup.broadinstitute.org/LDSCORE/
wget -P ~/uc-cross-atlas-data/ldsc \
  https://alkesgroup.broadinstitute.org/LDSCORE/eur_w_ld_chr.tar.bz2 \
  https://alkesgroup.broadinstitute.org/LDSCORE/w_hm3.snplist.bz2
cd ~/uc-cross-atlas-data/ldsc
tar -xjf eur_w_ld_chr.tar.bz2 && bunzip2 w_hm3.snplist.bz2

# Sanity check
python ~/ldsc/munge_sumstats.py --help | head
python ~/ldsc/ldsc.py --help | head
```

The wrappers in this module take `--ldsc-dir` and `--ld-dir` flags
that point at `~/ldsc` and `~/uc-cross-atlas-data/ldsc` respectively.

## Usage (after setup)

```bash
# Munge each GWAS to LDSC's .sumstats.gz format
python munge_for_ldsc.py \
    --in-snp-loc ../../data/gwas/uc_delange.snp.loc \
    --in-pval    ../../data/gwas/uc_delange.pval \
    --out        ../../results/ldsc/uc_delange.sumstats.gz \
    --trait UC --ldsc-dir ~/ldsc \
    --snplist ~/uc-cross-atlas-data/ldsc/w_hm3.snplist

# Run h2 + extract intercept / ratio
python run_ldsc_h2.py \
    --sumstats ../../results/ldsc/uc_delange.sumstats.gz \
    --ld-dir   ~/uc-cross-atlas-data/ldsc/eur_w_ld_chr \
    --out      ../../results/ldsc/uc_delange \
    --ldsc-dir ~/ldsc
```

`run_ldsc_h2.py` writes `uc_delange.log` (LDSC's full output) and
`uc_delange.intercept.tsv` (single-row: `gwas, h2, h2_se, intercept,
intercept_se, ratio, ratio_se, status`). The TSV is the F10 gate
deliverable.

Repeat for `yengo_height` (positive control), `scz_trubetskoy`
(negative control), and — after the (27)(c) ancestry-LD decision —
Liu.

## What lives in this module

- `README.md` — this file.
- `environment_py27.yml` — Python 2.7 conda env + LDSC's pinned
  dependencies (numpy 1.16, pandas 0.20, etc.).
- `munge_for_ldsc.py` — wrapper that joins `prepare_gwas.py` outputs
  back to a single sumstats file and calls LDSC's `munge_sumstats.py`.
- `run_ldsc_h2.py` — wrapper that calls `ldsc.py --h2` and parses
  the `.log` for the intercept + ratio.
- `parse_log.py` — LDSC `.log` parser (used by `run_ldsc_h2.py` and
  importable for tests).
- `test_munge_for_ldsc.py` — argument-construction tests.
- `test_run_ldsc_h2.py` — argument-construction + `.log` parsing
  tests.
- `test_parse_log.py` — `.log` parser unit tests against fixtures.

Tests run in the **`uc-cross-atlas` env** (the wrappers themselves
are Python 3); they do not require LDSC's Python 2.7 env to pass.
Execution against real GWAS happens in the Python 2.7 LDSC env.

## References

- Bulik-Sullivan et al. 2015 *Nat Genet* 47:291 (LDSC intercept).
- LDSC repo: `https://github.com/bulik/ldsc`.
- OPEN_FLAGS F10 (this module's gate).
- DECISIONS 23(a), 25, 27(a), 29(a) — N_eff discipline that feeds
  LDSC's `N` column.

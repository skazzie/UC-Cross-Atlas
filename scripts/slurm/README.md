# SLURM templates

One template per pipeline step. All scripts source `../config.sh`, so
edit that file (or export `UCC_*` env vars before submit) to set the
partition, account, and email.

## SBATCH defaults

All templates ship with the Hummingbird-verified defaults:

```
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --mail-user=amoli@ucsc.edu
```

To override at submit time (e.g., if a co-author runs the same scripts
under a different account):

```bash
sbatch -p $UCC_PARTITION --account=$UCC_ACCOUNT --mail-user=$UCC_EMAIL \
       --export=ALL,GWAS=delange \
       scripts/slurm/01_magma.slurm
```

## Catalog

| File | Purpose | Parameters | Notes |
|------|---------|------------|-------|
| `01_magma.slurm` | MAGMA pipeline for one GWAS | `GWAS=delange\|liu\|scz` | ~30 min, 16 GB. Submit 3 times. |
| `03_scdrs_compute.slurm` | scDRS regime-1 for one (atlas, GWAS) | `ATLAS=...`, `GWAS=...`, optional `GS_SUFFIX=_with_mhc`, `SEED=42` | ~30–45 min, 24 GB. Pan-GI: bump `--mem=48G`, `--time=06:00:00`, use highmem partition. |
| `04_seismic.slurm` | seismicGWAS regime-1 for one (atlas, GWAS) | `ATLAS=...`, `GWAS=...`, optional `PERMUTE=1` | ~30 min base; +30 min per (atlas,GWAS,tier) when `PERMUTE=1`. |
| `test_retest_array.slurm` | Test-retest seeds 1,2,3 for one (atlas, method) under de Lange | `ATLAS=...`, `METHOD=scdrs\|seismic`; `--array=1-3` baked in | 9 scDRS jobs + 9 seismic jobs total across the trio. |
| `donor_loo_array.slurm` | Donor-LOO jackknife for one (atlas, method) under de Lange | `ATLAS=...`, `METHOD=...`; set `--array=0-N` to match donor count | `garrido_trigo` has a built-in ≥5-cases-per-group guard. |

## M1 → M3 submission order

```bash
# 0. One-time setup (login node)
bash scripts/setup_env.sh
bash scripts/download_refs.sh
# (manually fetch atlases + Liu + SCZ per download_refs.sh notes)

# 1. M1: MAGMA for the three GWAS
sbatch --export=ALL,GWAS=delange scripts/slurm/01_magma.slurm
sbatch --export=ALL,GWAS=liu     scripts/slurm/01_magma.slurm
sbatch --export=ALL,GWAS=scz     scripts/slurm/01_magma.slurm

# 2. M3 regime-1 scDRS — 6 baseline runs (3 atlases x 2 GWAS)
for atlas in smillie garrido_trigo mennillo; do
  for gwas in delange liu; do
    sbatch --export=ALL,ATLAS=$atlas,GWAS=$gwas \
           scripts/slurm/03_scdrs_compute.slurm
  done
done

# 3. M3 regime-1 seismicGWAS — 6 baseline runs
for atlas in smillie garrido_trigo mennillo; do
  for gwas in delange liu; do
    sbatch --export=ALL,ATLAS=$atlas,GWAS=$gwas \
           scripts/slurm/04_seismic.slurm
  done
done

# 4. M3 sanity controls
sbatch --export=ALL,ATLAS=smillie,GWAS=scz scripts/slurm/03_scdrs_compute.slurm   # negative
sbatch --export=ALL,ATLAS=smillie,GWAS=delange,GS_SUFFIX=_with_mhc \
       scripts/slurm/03_scdrs_compute.slurm                                       # MHC sensitivity
# Tabula Muris x height (positive control) — separate atlas, write a small
# per-step submission once data is in place.

# 5. M3 test-retest — 3 atlases x 2 methods (de Lange only)
for atlas in smillie garrido_trigo mennillo; do
  for method in scdrs seismic; do
    sbatch --export=ALL,ATLAS=$atlas,METHOD=$method \
           scripts/slurm/test_retest_array.slurm
  done
done

# 6. M4 donor-LOO — broad tier, de Lange, both methods
sbatch --array=0-29 --export=ALL,ATLAS=smillie,METHOD=scdrs \
       scripts/slurm/donor_loo_array.slurm
sbatch --array=0-11 --export=ALL,ATLAS=garrido_trigo,METHOD=scdrs \
       scripts/slurm/donor_loo_array.slurm
sbatch --array=0-N  --export=ALL,ATLAS=mennillo,METHOD=scdrs \
       scripts/slurm/donor_loo_array.slurm
# (and the seismic equivalents)
```

Replace `0-N` with the actual Mennillo donor count after pre-treatment
subsetting (verify ≥8 per DECISIONS.md).

## Monitoring

```bash
squeue -u $USER                          # pending + running
squeue -u $USER -t R                     # running only
sacct -u $USER --starttime=$(date -d '1 day ago' +%Y-%m-%d) \
      --format=JobID,JobName,State,Elapsed,MaxRSS,ExitCode

scontrol show job <jobid>                # detail of one job
seff <jobid>                             # post-hoc efficiency report
```

Logs land in the directory you submitted from (because of `--output=%x_%j.out`).
Move them or change the `--output=` path if you'd rather they go to
`$UCC_LOGS`.

## Things that aren't templated yet

- Pan-GI scDRS — same as `03_scdrs_compute.slurm` but `--mem=48G`,
  longer time, highmem partition. Submit with overrides rather than a
  separate file.
- Brown's method (`code/07_regime2_meta/`) — pure post-processing, fast
  enough to run interactively or with a small SLURM job; M5 work.
- HCA Gut + Pan-GI broad-comparator concordance — same `03_scdrs_compute`
  template, just different `ATLAS=` values; M5 work.
- Tabula Muris × height (positive control) — `03_scdrs_compute` with
  `ATLAS=tabula_muris,GWAS=yengo_height` once those files exist.

## If a job fails

1. Check the `.err` file. Most failures are `ERROR: missing input ...`
   (bad path) or OOM (re-submit with `--mem=Xg` higher).
2. For OOM kills, `seff <jobid>` shows `MaxRSS`. Set `--mem` to roughly
   1.3× that.
3. For `module: command not found` or conda failures, edit
   `setup_env.sh` and the `module load` lines in each `.slurm` to match
   what `module avail` shows on Hummingbird.

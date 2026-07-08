# Chains scDRS compute-score + perform-downstream for the two control
# GWAS on Garrido. Run in background; logs to AppData\Local\Temp.
#
# Sequence (MHC-excluded only — matches locked v1 primary, the
# apples-to-apples comparison against UC):
#   1. SCZ compute-score      ~6 min
#   2. SCZ pdownstream broad  ~13 min
#   3. Height compute-score   ~6 min
#   4. Height pdownstream broad ~13 min
# Total wall ~38 min. MHC-included variants of controls deferred.

Set-Location 'C:\Users\muska\UC-Cross-Atlas'
$ErrorActionPreference = 'Continue'

function StepLog {
  param([string]$msg)
  $ts = (Get-Date).ToString("HH:mm:ss")
  Write-Output "[$ts] $msg"
}

# ---- SCZ compute-score is launched OUTSIDE this script (already running).
# Wait for its output before pdownstream.

$scz_full = 'results\scdrs\garrido_scz\SCZ.full_score.gz'
StepLog "Waiting for $scz_full ..."
while (-not (Test-Path $scz_full)) { Start-Sleep -Seconds 30 }
StepLog "SCZ full_score landed; running pdownstream"

& py.exe scdrs_cli.py perform-downstream `
  'data\atlases\garrido.h5ad' `
  'results\scdrs\garrido_scz\@.full_score.gz' `
  'results\scdrs\garrido_scz\' `
  --group_analysis=cell_type_broad `
  --flag_filter_data=True --flag_raw_count=False 2>&1 | `
  Tee-Object 'C:\Users\muska\AppData\Local\Temp\scdrs_scz_pdown.log'
StepLog "SCZ pdownstream done (exit=$LASTEXITCODE)"

# ---- Height compute-score
StepLog "Launching Height compute-score"
& py.exe scdrs_cli.py compute-score `
  --h5ad_file data/atlases/garrido.h5ad --h5ad_species human `
  --gs_file results/magma/yengo_height.gs --gs_species human `
  --cov_file data/atlases/garrido_covariates.tsv `
  --flag_filter_data True --flag_raw_count False --n_ctrl 1000 `
  --out_folder results/scdrs/garrido_height/ 2>&1 | `
  Tee-Object 'C:\Users\muska\AppData\Local\Temp\scdrs_height_compute.log'
StepLog "Height compute-score done (exit=$LASTEXITCODE)"

# ---- Height pdownstream
StepLog "Launching Height pdownstream"
& py.exe scdrs_cli.py perform-downstream `
  'data\atlases\garrido.h5ad' `
  'results\scdrs\garrido_height\@.full_score.gz' `
  'results\scdrs\garrido_height\' `
  --group_analysis=cell_type_broad `
  --flag_filter_data=True --flag_raw_count=False 2>&1 | `
  Tee-Object 'C:\Users\muska\AppData\Local\Temp\scdrs_height_pdown.log'
StepLog "Height pdownstream done (exit=$LASTEXITCODE)"

StepLog "ALL CONTROLS DONE"

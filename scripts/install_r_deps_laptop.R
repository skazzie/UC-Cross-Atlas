# Laptop-side R dep installer for the uc-cross-atlas-r env.
#
# Run inside the env via:
#   conda run -n uc-cross-atlas-r Rscript scripts/install_r_deps_laptop.R
#
# Differs from scripts/install_r_deps.R in one way: that script assumes
# the bioconductor stack was conda-installed (it stops loudly if missing).
# This script installs the bioconductor packages itself via BiocManager,
# because the laptop env (scripts/environment-r-laptop.yml) deliberately
# does not pull the bioconductor stack from conda — bioconda's Windows
# builds either don't exist or conflict.

repos <- "https://cloud.r-project.org"

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = repos)
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = repos)
}

# Bioconductor packages seismicGWAS + downstream code (07_regime2_meta,
# 04_seismic h5ad loader) require.
bioc_pkgs <- c(
  "SingleCellExperiment",
  "EmpiricalBrownsMethod",
  "zellkonverter"
)
need <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(need)) {
  cat(sprintf("[install_r_deps_laptop] installing via BiocManager: %s\n",
              paste(need, collapse = ", ")))
  BiocManager::install(need, update = FALSE, ask = FALSE)
}

# CRAN extras needed by run_seismic.R / install_r_deps.R.
cran_pkgs <- c("ontologyIndex")
need <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(need)) {
  install.packages(need, repos = repos)
}

# seismicGWAS + scOntoMatch from GitHub.
gh_pkgs <- list(
  seismicGWAS  = "ylaboratory/seismicGWAS",
  scOntoMatch  = "Papatheodorou-Group/scOntoMatch"
)
for (pkg in names(gh_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[install_r_deps_laptop] installing %s from GitHub\n", pkg))
    remotes::install_github(gh_pkgs[[pkg]], upgrade = "never")
  } else {
    cat(sprintf("[install_r_deps_laptop] %s already installed\n", pkg))
  }
}

cat("\n[install_r_deps_laptop] done\n")
for (p in c("SingleCellExperiment", "EmpiricalBrownsMethod", "zellkonverter",
            "ontologyIndex", "seismicGWAS", "scOntoMatch")) {
  cat(sprintf("  %-22s %s\n", p, requireNamespace(p, quietly = TRUE)))
}

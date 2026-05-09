# R packages not on conda (seismicGWAS, scOntoMatch from GitHub).
# Run AFTER setup_env.sh has activated the conda env, so the conda R is
# what `Rscript` resolves to.

repos <- "https://cloud.r-project.org"

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = repos)
}

# Optional sanity: fail loudly if conda's seismicGWAS prerequisites are missing.
required <- c(
  "SingleCellExperiment", "EmpiricalBrownsMethod",
  "ontologyIndex", "tidyverse", "remotes"
)
missing <- setdiff(required, rownames(installed.packages()))
if (length(missing)) {
  stop(
    "Missing prerequisite R packages (should have been installed by conda): ",
    paste(missing, collapse = ", "),
    "\nCheck scripts/environment.yml and re-run setup_env.sh."
  )
}

gh_pkgs <- list(
  seismicGWAS = "ylaboratory/seismicGWAS",
  scOntoMatch = "Papatheodorou-Group/scOntoMatch"
)

for (pkg in names(gh_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[install_r_deps] installing %s from GitHub\n", pkg))
    remotes::install_github(gh_pkgs[[pkg]], upgrade = "never")
  } else {
    cat(sprintf("[install_r_deps] %s already installed\n", pkg))
  }
}

cat("\n[install_r_deps] done\n")

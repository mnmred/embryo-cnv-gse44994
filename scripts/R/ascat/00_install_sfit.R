# ASCAT trial, step 0 - install the one R package that conda does not have.
#
# CRMA v2's first step (allelic crosstalk calibration) fits a cone to the A/B intensities
# with the package sfit, by the author of aroma.affymetrix. sfit is neither on conda-forge,
# bioconda nor CRAN, so it is built from the author's source, pinned to release 0.3.2.
#
# Run once (ascat_env active):   Rscript scripts/R/ascat/00_install_sfit.R

SFIT <- "https://github.com/HenrikBengtsson/sfit/archive/refs/tags/0.3.2.tar.gz"

if (!requireNamespace("sfit", quietly = TRUE)) {
  install.packages(SFIT, repos = NULL, type = "source")
}
message("sfit ", as.character(packageVersion("sfit")))

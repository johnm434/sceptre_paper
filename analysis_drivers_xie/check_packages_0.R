# This script verifies that all packages required for the Xie analysis are available.
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/utilities/verify_all_packages_available.R"))

packages <- c("R.matlab", "R.utils", "fst", "sn", "MASS", "VGAM", "tidyverse", "bigstatsr", "openxlsx", "rhdf5", "ravel", "sceptre")
locs <- c(rep("CRAN", 9), "Bioc", rep("github", 2))
github_repo <- c(rep(NA, 10), "Timothy-Barry/ravel", "Timothy-Barry/sceptre_paper")
github_repo_subdir <- c(rep(NA, 11), "sceptre")
df <- data.frame(package = packages, loc = locs, github_repo = github_repo, github_repo_subdir = github_repo_subdir)

verify_all_packages_available(df)

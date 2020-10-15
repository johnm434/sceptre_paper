# A function to ensure that all packages are available. If not, the function downloads and installs missing packages.
# df is a data frame with four columns: (i) package name, (ii) package location (one of "CRAN", "Bioc", and "github"), (iii) github repository (github only), and (iv) github sub-repository (github only)

verify_all_packages_available <- function(df) {
  my_packs <- rownames(installed.packages())
  for (i in 1:nrow(df)) {
    curr_package <- df[i, "package"]
    cat(paste("Checking", curr_package, "\n"))
    if(!(curr_package %in% my_packs)) {
      curr_loc <- df[i, "loc"]
      if (curr_loc == "CRAN")  {
        install.packages(curr_package)
      } else if (curr_loc == "Bioc") {
        if (!("BiocManager" %in% my_packs)) install.packages("BiocManager")
        library(BiocManager)
        BiocManager::install(curr_package)
      } else if (curr_loc == "github") {
        git_dir <- df[i, "github_repo"]
        git_subdir <- df[i, "github_repo_subdir"]
        if (!("devtools" %in% my_packs)) install.packages("devtools")
        library(devtools)
        if (is.na(git_subdir)) install_github(repo = git_dir) else install_github(repo = git_dir, subdir = git_subdir)
      }
    }
  }
}

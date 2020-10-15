# A function to ensure that all packages are available. If not, the function downloads and installs missing packages.
# df is a data frame with four columns: (i) package name, (ii) package location (one of "CRAN", "Bioc", and "github"), (iii) github repository (github only), and (iv) github sub-repository (github only)

verify_all_packages_available <- function(df) {
  my_packs <- rownames(installed.packages())
  for (i in 1:nrow(df)) {
    curr_package <- as.character(df[[i, "package"]])
    cat(paste("Checking", curr_package, "\n"))
    if(!(curr_package %in% my_packs)) {
      curr_loc <- as.character(df[[i, "loc"]])
      if (curr_loc == "CRAN")  {
        install.packages(curr_package, repos = "https://cloud.r-project.org")
      } else if (curr_loc == "Bioc") {
        if (!("BiocManager" %in% my_packs)) install.packages("BiocManager")
        library(BiocManager)
        BiocManager::install(curr_package)
      } else if (curr_loc == "github") {
        git_dir <- as.character(df[[i, "github_repo"]])
        git_subdir <- as.character(df[[i, "github_repo_subdir"]])
        if (!("devtools" %in% my_packs)) install.packages("devtools")
        library(devtools)
        if (is.na(git_subdir)) install_github(repo = git_dir) else install_github(repo = git_dir, subdir = git_subdir)
      }
    }
  }
}

# A small utility function to create strings representing all parent directories of s
create_parent_directories <- function(s) {
  dirs <- unlist(str_split(string = s, pattern = "/"))
  out <- map(.x = 1:length(dirs), .f = function(i) paste0(dirs[1:i], collapse = "/"))
  return(unlist(out))
}

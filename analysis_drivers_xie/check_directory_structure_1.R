args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers_xie/paths_to_dirs.R"))
packages <- c("purrr", "stringr")
for (package in packages) suppressPackageStartupMessages(library(package, character.only = TRUE))

# Hardcode the directories to create.
sub_dirs <- c(create_parent_directories("data/xie/raw"), create_parent_directories("data/xie/precomp/gRNA"), "data/xie/precomp/gene", "data/xie/processed",
              create_parent_directories("results/xie"), "figures", create_parent_directories("logs/xie")) %>% unique()

dirs_to_create <- paste0(offsite_dir, "/", sub_dirs)
for (directory in dirs_to_create) {
  if (!dir.exists(directory)) {
    dir.create(directory)
  }
}

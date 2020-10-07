args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
packages <- c("purrr", "stringr")
for (package in packages) suppressPackageStartupMessages(library(package, character.only = TRUE))

# Hardcode the directories to create.
sub_dirs <- c(create_parent_directories("data/gasperini/raw"), create_parent_directories("data/gasperini/precomp/gRNA"), "data/gasperini/precomp/gene", "data/gasperini/processed",
  create_parent_directories("data/functional"), "data/functional/HIC", "data/functional/ChIP-seq", "data/functional/GeneHancer",
  create_parent_directories("results/gasperini"), "results/simulations",
  "figures", create_parent_directories("logs/gasperini")) %>% unique()

dirs_to_create <- paste0(offsite_dir, "/", sub_dirs)
for (directory in dirs_to_create) {
  if (!dir.exists(directory)) {
    dir.create(directory)
  }
}

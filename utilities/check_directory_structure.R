library(stringr, quietly = TRUE)
library(purrr, quietly = TRUE)
# take as an argument a file path to the "offsite" (data, results, logs, figures) directory.
args <- commandArgs(trailingOnly = TRUE)
offsite_directory <- args[1]

# This function takes as an argument an inner directory, and returns all 
create_parent_directories <- function(s) {
  s_split <- str_split(string = s, pattern = "/")[[1]]
  unlist(map(1:length(s_split), function(i) paste0(s_split[1:i], collapse = "/")))
}

# Hardcode the directories to create.
sub_dirs <- c(create_parent_directories("data/gasperini/raw"), "data/gasperini/precomp", "data/gasperini/processed",
  create_parent_directories("data/Xie/raw"), "data/Xie/processed", "data/Xie/precomp",
  create_parent_directories("data/functional"), "data/functional/HIC", "data/functional/ChIP-seq", "data/functional/GeneHancer", 
  create_parent_directories("results/gasperini"), "results/Xie", "results/simulations",
  "figures", "logs") %>% unique()

dirs_to_create <- paste0(offsite_directory, "/", sub_dirs)

for (directory in dirs_to_create) {
  if (!dir.exists(directory)) {
    dir.create(directory)
  }
}

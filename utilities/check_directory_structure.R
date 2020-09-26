library(stringr, quietly = TRUE)
args <- commandArgs(trailingOnly = TRUE)
data_res_directory <- args[1]

data_res_sub_directories <- c("data", "data/raw", "data/raw/CRISPR", "data/raw/ChIP-seq",
                    "data/raw/HIC", "data/raw/GeneHancer", "data/processed", "precomp",
                    "results", "results/pvalues", "results/resampled_zvalues",
                    "figures")

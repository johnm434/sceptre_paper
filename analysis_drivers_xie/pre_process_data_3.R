# Pre-process data
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers_xie/paths_to_dirs.R"))
suppressPackageStartupMessages(library(rhdf5))

# First, we pre-process the expression matrices; this involves loading each expression matrix, extracting the UMI counts, and storing in a file-backed matrix.
h5_obj <- H5Fopen(paste0(raw_data_dir, "/GSM3722727_K562-dCas9-KRAB_5K-sgRNAs_Batch-4_1_filtered_gene_bc_matrices_h5.h5"))


args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
# Define the code and offsite dirs
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
param_file <- if(is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/analysis_drivers_xie/sceptre_function_args.R" else args[2]
# source the function arguments
source(param_file)

pod_id <- as.integer(args[3])
run_gene_precomputation_at_scale(pod_id = pod_id, gene_precomp_dir = gene_precomp_dir, cell_gene_expression_matrix = cell_gene_expression_matrix, ordered_gene_ids = ordered_gene_ids, covariate_matrix = covariate_matrix, cell_subset = cell_subset, gene_sizes = gene_sizes, log_dir = log_dir)

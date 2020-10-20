args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
suppressPackageStartupMessages(library(katsevich2020))
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
param_file <- if(is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/analysis_drivers_xie/sceptre_function_args.R" else args[2]
source(param_file)

pod_id <- as.integer(args[3])
map_reduce_at_scale(f = "run_NB_model", pod_id = pod_id, dictionary = paste0(offsite_dir, "/results/xie/dictionary.fst") %>% read.fst, results_dir = results_dir_negbin, cell_gene_expression_matrix = cell_gene_expression_matrix, ordered_gene_ids = ordered_gene_ids, gRNA_indicator_matrix_fp = gRNA_indicator_matrix_fp, covariate_matrix = covariate_matrix, cell_subset = cell_subset, log_dir = log_dir)

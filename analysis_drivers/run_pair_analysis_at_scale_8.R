# Run pair analysis at scale
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers/file_paths_to_dirs.R"))
source(paste0(code_dir, "/analysis_drivers/sceptre_function_args.R"))

pod_id <- as.integer(args[3]) # pod input
run_gRNA_gene_pair_analysis_at_scale(pod_id, gene_precomp_dir, gRNA_precomp_dir, results_dir, cell_gene_expression_matrix, ordered_gene_ids, gRNA_indicator_matrix_fp, covariate_matrix, cell_subset, log_dir, seed = 1234)

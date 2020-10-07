# run gene precomputation at scale
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers/file_paths_to_dirs.R"))
source(paste0(code_dir, "/analysis_drivers/sceptre_function_args.R"))

pod_id <- as.integer(args[3]) # pod_id input
run_gene_precomputation_at_scale(pod_id = pod_id, gene_precomp_dir = gene_precomp_dir, cell_gene_expression_matrix = cell_gene_expression_matrix, ordered_gene_ids = ordered_gene_ids, covariate_matrix = covariate_matrix, cell_subset = cell_subset, gene_sizes = select_sizes, log_dir = log_dir)

# Run gRNA precomputation at scale
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers_gasp/file_paths_to_dirs.R"))
source(paste0(code_dir, "/analysis_drivers_gasp/sceptre_function_args.R"))

pod_id <- as.integer(args[3]) # pod_id input
run_gRNA_precomputation_at_scale(pod_id = pod_id, gRNA_precomp_dir = gRNA_precomp_dir, gRNA_indicator_matrix_fp = gRNA_indicator_matrix_fp, covariate_matrix = covariate_matrix, cell_subset = cell_subset, log_dir = log_dir)

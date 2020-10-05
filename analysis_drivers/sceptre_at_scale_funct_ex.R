# Run pair analysis at scale
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers/file_paths_to_dirs.R"))
source(paste0(code_dir, "/analysis_drivers/sceptre_function_args.R"))
gRNA_gene_pairs <- read.fst(paste0(processed_dir, "/gene_gRNA_pairs_to_study.fst"))
gRNA_gene_pairs <- gRNA_gene_pairs %>% slice(sample(1:nrow(gRNA_gene_pairs), 10))

res <- run_sceptre_at_scale(gRNA_gene_pairs = gRNA_gene_pairs,
                     gene_precomp_dir = gene_precomp_dir,
                     gRNA_precomp_dir = gRNA_precomp_dir,
                     results_dir = results_dir,
                     cell_gene_expression_matrix = cell_gene_expression_matrix,
                     ordered_gene_ids = ordered_gene_ids,
                     gRNA_indicator_matrix_fp = gRNA_indicator_matrix_fp,
                     covariate_matrix = covariate_matrix,
                     cell_subset = cell_subset,
                     pod_sizes = c(gene = 5, gRNA = 5, pair = 5),
                     seed = 1234,
                     log_dir = log_dir,
                     multi_processor = TRUE)

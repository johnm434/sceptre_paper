args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers_xie/paths_to_dirs.R"))

exp_mat_t <- readRDS(paste0(processed_dir, "/exp_mat_t_metadata.rds")) %>% load_fbm()
# n_genes_per_cell <- big_apply(exp_mat_t, function(X, ind) {colSums(X[,ind] > 0)}) %>% unlist()
n_umis_per_cell <- big_apply(exp_mat_t, function(X, ind) {colSums(X[,ind])}) %>% unlist()

covariate_matrix <- read.fst(paste0(processed_dir, "/cell_covariate_matrix.fst"))
covariate_matrix <- covariate_matrix %>% mutate(n_umis = n_umis_per_cell)
covariate_model_matrix <- summarize(covariate_matrix, batch = paste0("batch_", batch) %>% factor, log_n_umis = log(n_umis))
write.fst(x = covariate_model_matrix, path = paste0(processed_dir, "/covariate_model_matrix.fst"))

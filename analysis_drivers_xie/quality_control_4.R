args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers_xie/paths_to_dirs.R"))

# First, save the model covariate matrix
exp_mat_t <- readRDS(paste0(processed_dir, "/exp_mat_t_metadata.rds")) %>% load_fbm()
# n_genes_per_cell <- big_apply(exp_mat_t, function(X, ind) {colSums(X[,ind] > 0)}) %>% unlist()
n_umis_per_cell <- big_apply(exp_mat_t, function(X, ind) {colSums(X[,ind])}) %>% unlist()
covariate_matrix <- read.fst(paste0(processed_dir, "/cell_covariate_matrix.fst"))
covariate_model_matrix <- covariate_matrix %>% mutate(n_umis = n_umis_per_cell) %>% summarize(batch = paste0("batch_", batch) %>% factor, log_n_umis = log(n_umis))
write.fst(x = covariate_model_matrix, path = paste0(processed_dir, "/covariate_model_matrix.fst"))

# Next, determine which genes to analyze
gRNA_indic_mat <- read.fst(paste0(processed_dir, "/gRNA_indicator_matrix.fst"))
gRNA_id <- colnames(gRNA_indic_mat)
gene_ids <- readRDS(paste0(processed_dir, "/ordered_genes.RDS"))
exp_mat <- readRDS(paste0(processed_dir, "/exp_mat_metadata.rds")) %>% load_fbm()
gene_expression_p <- big_apply(exp_mat, function(X, ind) colMeans(X[,ind] >= 1)) %>% unlist()
highly_expressed_genes <- gene_ids[which(gene_expression_p >= 0.05)]
gRNA_gene_pairs <- tibble(gene_id = highly_expressed_genes, gRNA_id = gRNA_id)
write.fst(gRNA_gene_pairs, paste0(processed_dir, "/gRNA_gene_pairs.fst"))
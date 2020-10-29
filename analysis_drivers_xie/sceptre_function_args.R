# sceptre function arguments; these arguments should be defined in terms of "offsite_dir"
# offsite_dir <- "/Volumes/tims_new_drive/research/sceptre_files"
small_example <- TRUE

processed_dir <- paste0(offsite_dir, "/data/xie/processed")
gene_precomp_dir <- paste0(offsite_dir, "/data/xie/precomp/gene")
gRNA_precomp_dir <- paste0(offsite_dir, "/data/xie/precomp/gRNA")
results_dir <- paste0(offsite_dir, "/results/xie/sceptre")
results_dir_negbin <- paste0(offsite_dir, "/results/xie/negative_binomial")
log_dir <- paste0(offsite_dir, "/logs/xie")
gRNA_gene_pairs <- read.fst(paste0(processed_dir, "/gRNA_gene_pairs.fst"))
covariate_matrix <- read.fst(paste0(processed_dir, "/covariate_model_matrix.fst"))
cell_gene_expression_matrix <- readRDS(paste0(processed_dir, "/exp_mat_metadata.rds")) %>% load_fbm
ordered_gene_ids <- readRDS(paste0(processed_dir, "/ordered_gene_ids.RDS"))
gRNA_indicator_matrix_fp <- paste0(processed_dir, "/gRNA_indicator_matrix.fst")
regularization_amount <- 3
cell_subset <- readRDS(paste0(processed_dir, "/cell_subsets.rds"))[["exploratory_cells"]]
seed <- 1234
B <- 500
pod_sizes <- c(gene = 200, gRNA = 1, pair = 200)
gene_sizes <- NULL

if (small_example) {
  pod_sizes <- c(gene = 10, gRNA = 1, pair = 10)
  gRNA_gene_pairs <- slice(gRNA_gene_pairs, 1:20)
}

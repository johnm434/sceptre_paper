# sceptre function arguments

covariate_matrix <- read.fst(paste0(processed_dir, "/cell_covariate_model_matrix.fst"))
cell_subset <- readRDS(paste0(processed_dir, "/cells_to_keep.rds"))
cell_gene_expression_matrix <- readRDS(paste0(processed_dir, "/expression_FBM_metadata.rds")) %>% load_fbm
ordered_gene_ids <- readRDS(paste0(processed_dir, "/ordered_gene_ids.RDS"))
gRNA_indicator_matrix_fp <- paste0(processed_dir, "/gRNA_indicators.fst")

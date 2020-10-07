# sceptre function arguments

covariate_matrix <- read.fst(paste0(processed_dir, "/cell_covariate_model_matrix.fst"))
cell_subset <- readRDS(paste0(processed_dir, "/cells_to_keep.rds"))
cell_gene_expression_matrix <- readRDS(paste0(processed_dir, "/expression_FBM_metadata.rds")) %>% load_fbm
ordered_gene_ids <- readRDS(paste0(processed_dir, "/ordered_gene_ids.RDS"))
gRNA_indicator_matrix_fp <- paste0(processed_dir, "/gRNA_indicators.fst")
gRNA_gene_pairs <- read.fst(paste0(processed_dir, "/gene_gRNA_pairs_to_study.fst"))
if (TRUE) gRNA_gene_pairs <- gRNA_gene_pairs %>% slice(c(12262, 185202, 323756, 305770, 492670, 596203, 553777, 205062, 400245, 320979, 288925, 506135, 607659, 222201, 332665, 519542, 614764, 457555, 358113, 481437, 166854, 304582, 141405, 360173, 64637))
all_dispersions <- read.fst(paste0(processed_dir, "/disp_table.fst"))
select_dispersions <- sapply(X = gRNA_gene_pairs$gene_id, FUN = function(id) {
  filter(all_dispersions, gene_id == as.character(id)) %>% pull(disp)
  })
select_sizes <- 1/select_dispersions
names(select_sizes) <- gRNA_gene_pairs$gene_id
rm(all_dispersions, select_dispersions)

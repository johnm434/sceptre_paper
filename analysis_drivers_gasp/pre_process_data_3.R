args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers/file_paths_to_dirs.R"))

#######################
# Gasperini CRISPR data
#######################

# 1. Create file-backed matrices for the .mtx expression matrix; both the gene-cell matrix and its transpose.
path_to_mtx <- paste0(raw_data_dir, "/GSE120861_at_scale_screen.exprs.mtx")
fbms <- create_fbm_from_mtx(path_to_mtx, processed_dir)
# save the .bk metadata for both of these matrices (nrow, ncell, data type, file location)
for (i in 1:length(fbms)) {
  file_name <- paste0(names(fbms)[[i]], ".rds")
  full_file_path <- paste0(processed_dir, "/",file_name)
  saveRDS(fbms[[i]], full_file_path)
}

# 2. Read monocole object
gc()
library(monocle, quietly = TRUE) # monocole required to load monocole CellDataSet object; we will extract only the cell-specific metadata
m_object <- readRDS(paste0(raw_data_dir, "/GSE120861_at_scale_screen.cds.rds"))

# 2a. save the coefficients of the mean-dispersion relationship
saveRDS(attr(m_object@dispFitInfo$blind$disp_func, "coefficients"), paste0(processed_dir, "/disp_coefficients.rds"))

# 2b. save the raw dispersions
write.fst(m_object@dispFitInfo$blind$disp_table, paste0(processed_dir, "/disp_table.fst"))

cell_metadata <- pData(m_object)
rm(m_object); gc()
covariates_cols <- 1:18
gRNA_cols <- 19:ncol(cell_metadata)

# 2c. Save the gRNA indicators
gRNA_indicators <- cell_metadata[,gRNA_cols]
write.fst(x = gRNA_indicators, compress = 0, path = paste0(processed_dir, "/gRNA_indicators.fst"))

# 2d. Save the covariates
cell_covariates <- cell_metadata[,covariates_cols]
write.fst(x = cell_covariates, path = paste0(processed_dir, "/cell_covariates_all.fst"))

# 3. Create and save the covariate matrix that we will use in the regressions.
cell_covariates <- read.fst(paste0(processed_dir, "/cell_covariates_all.fst"))

# Identify the first 5 variables and apply appropriate transformations
cell_covariates_model <- summarize(cell_covariates, p_mito = cell_covariates$percent.mito, prep_batch = cell_covariates$prep_batch, lg_total_umis = log(cell_covariates$total_umis), lg_guide_count = log(cell_covariates$guide_count))

# Determine the number of expressed genes in each cell
exp_mat_t <- readRDS(paste0(processed_dir, "/expression_FBM_t_metadata.rds")) %>% load_fbm
n_genes <- big_apply(exp_mat_t, function(X, ind) {colSums(X[,ind] > 0)}) %>% unlist()

# Add this feature to the cell_covariates_model data frame
cell_covariates_model <- mutate(cell_covariates_model, lg_n_genes = log(n_genes))
write.fst(cell_covariates_model, paste0(processed_dir, "/cell_covariate_model_matrix.fst"))

# 4. Extract the target_site-gRNA pairs to analyze.
all_deg_results <- suppressWarnings(read_tsv(paste0(raw_data_dir, "/GSE120861_all_deg_results.at_scale.txt"), col_types = "cddddddccccciiciiccl"))
pairs_to_analyze <- all_deg_results %>% rename(gene_id = ENSG, gRNA_id = gRNA_group) %>% select(gene_id, gRNA_id) %>% mutate(gene_id = factor(gene_id), gRNA_id = factor(gRNA_id)) %>% arrange()
write.fst(pairs_to_analyze, paste0(processed_dir, "/gene_gRNA_pairs_to_study.fst"))

# 5. Finally, save the (ordered) gene names
gene_ids <- read_tsv(paste0(raw_data_dir, "/GSE120861_at_scale_screen.genes.txt"),
                    col_names = FALSE, col_types = "c") %>% pull()
saveRDS(gene_ids, paste0(processed_dir, "/ordered_gene_ids.RDS"))

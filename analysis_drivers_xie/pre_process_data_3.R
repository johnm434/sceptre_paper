# Pre-process data
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers_xie/paths_to_dirs.R"))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(bigstatsr))
suppressPackageStartupMessages(library(stringr))

###################
# Expression matrix
###################

# First, we determine the number of cells and number of genes across all the batches
h5_files <- paste0(raw_data_dir, "/", grep(pattern = '*.h5', list.files(raw_data_dir), value = TRUE))
dims_across_h5s <- sapply(h5_files, function(h5_file) {
  h5_handle <- H5Fopen(h5_file)
  dim <- h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/shape"
  out <- as.integer(dim)
  H5Fclose(h5_handle)
  return(out)
}) %>% t()
row.names(dims_across_h5s) <- NULL
all(dims_across_h5s[,1] == dims_across_h5s[1,1]) # Verify n genes consistent across files
n_cells_total <- sum(dims_across_h5s[,2])
n_genes_total <- dims_across_h5s[1,1]

# Select a subset of genes for which to store the data
n_genes_in_use <- 1000
genes_to_use_idxs <- sample(1:n_genes_total, n_genes_in_use) %>% sort()

# Next, we create a file-backed matrix to store the transpose of the expression matrix
exp_mat_t <- FBM(nrow = n_genes_in_use, ncol = n_cells_total, type = "unsigned short", init = 0, backingfile = paste0(processed_dir, "/expression_matrix_t"), create_bk = TRUE)

# We iterate through the hd5 files and write the column chunks piece-by-piece to the FBM.
n_cells_processed <- 0
batch <- integer(0)
for (h5_file in h5_files) {
  print(paste("Working on", h5_file))
  h5_handle <- H5Fopen(h5_file)
  dat <- h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/data"
  indices <- h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/indices"
  ind_ptr <- h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/indptr"
  shape <-  h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/shape"
  n_rows <- shape[1]
  n_cols <- shape[2]
  to_write <- extract_column_from_csc(column_no = 1:n_cols, dat = dat, indices = indices, ind_ptr = ind_ptr, n_rows = n_rows, row_idxs = genes_to_use_idxs)
  exp_mat_t[1:n_genes_in_use, (n_cells_processed + 1):(n_cells_processed + n_cols)] <- to_write
  H5Fclose(h5_handle)
  n_cells_processed <- n_cells_processed + n_cols
  curr_batch <- str_extract(string = h5_file, pattern = "Batch-[0-9]") %>% str_extract(pattern = "[0-9]") %>% as.integer()
  batch <- c(batch, rep(curr_batch, n_cols))
}

exp_mat <- big_transpose(exp_mat_t, backingfile = paste0(processed_dir, "/expression_matrix"))
cell_covariate_matrix <- tibble(batch = batch)

##############
# Bulk RNA-seq
##############
read_tsv(paste0(raw_data_dir, "/GSE129825_Libraries.FeatureCounts.ARL15_enhancer.txt"))

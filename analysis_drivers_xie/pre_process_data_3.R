# Pre-process data
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers_xie/paths_to_dirs.R"))
packs <- c("rhdf5", "stringr", "openxlsx", "ravel")
for (pack in pakcs) suppressPackageStartupMessages(library(pack, character.only = TRUE))

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

# We will use only the protein-coding genes in this analysis; load them.
gene_df <- read.xlsx(xlsxFile = paste0(raw_data_dir, "/Genes.xlsx"), sheet = 1)
all_protein_coding_genes <- gene_df$Gene_Symbol
rm(gene_df)
h5_file <- h5_files[1]
h5_handle <- H5Fopen(h5_file)
all_sequenced_genes <- h5_handle$"/refgenome_hg38_CROP-Guide-MS2-2.1.0/gene_names"
genes_in_use <- all_sequenced_genes[which(all_sequenced_genes %in% all_protein_coding_genes)]
n_genes_in_use <- length(genes_in_use)
H5Fclose(h5_handle)

# Next, we create a file-backed matrix to store the transpose of the expression matrix
exp_mat_t <- FBM(nrow = n_genes_in_use, ncol = n_cells_total, type = "unsigned short", init = 0, backingfile = paste0(processed_dir, "/expression_matrix_t"), create_bk = TRUE)
exp_mat_t_metadata <- list(nrow = n_genes_in_use, ncol = n_cells_total, type = "unsigned short", backingfile = paste0(processed_dir, "/expression_matrix_t"))
saveRDS(object = exp_mat_t_metadata, file = paste0(processed_dir, "/exp_mat_t_metadata.rds"))

exp_mat <- FBM(nrow = n_cells_total, ncol = n_genes_in_use, type = "unsigned short", init = 0, backingfile = paste0(processed_dir, "/expression_matrix"), create_bk = TRUE)
exp_mat_metadata <- list(nrow = n_cells_total, ncol = n_genes_in_use, type = "unsigned short", backingfile = paste0(processed_dir, "/expression_matrix"))
saveRDS(object = exp_mat_metadata, file = paste0(processed_dir, "/exp_mat_metadata.rds"))

# We iterate through the hd5 files and write the column chunks piece-by-piece to the FBM.
n_cells_processed <- 0
batch <- integer(0)
ordered_cell_barcodes <- character(0)
for (h5_file in h5_files) {
  print(paste("Working on", h5_file))
  h5_handle <- H5Fopen(h5_file)
  dat <- h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/data"
  indices <- h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/indices"
  ind_ptr <- h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/indptr"
  shape <-  h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/shape"
  n_cols <- shape[2]
  all_gene_names <- h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/gene_names"
  cell_barcodes <- h5_handle$"refgenome_hg38_CROP-Guide-MS2-2.1.0/barcodes"
  row_subset <- match(x = genes_in_use, table = all_gene_names)
  # Ensure that the maximum index pointer is less than the maximum integer value in R
  if (max(ind_ptr) >= 2147483647) stop("Index pointer exceeds R max integer value.")
  to_write <- extract_column_from_csc_matrix(col_nos = 1:n_cols, csc_mat = list(p = ind_ptr, i = indices, x = dat, dim = shape), row_subset = row_subset)
  exp_mat_t[1:n_genes_in_use, (n_cells_processed + 1):(n_cells_processed + n_cols)] <- to_write
  exp_mat[(n_cells_processed + 1):(n_cells_processed + n_cols), 1:n_genes_in_use] <- t(to_write)
  H5Fclose(h5_handle)
  n_cells_processed <- n_cells_processed + n_cols
  curr_batch <- str_extract(string = h5_file, pattern = "Batch-[0-9]") %>% str_extract(pattern = "[0-9]") %>% as.integer()
  batch <- c(batch, rep(curr_batch, n_cols))
  ordered_cell_barcodes <- c(ordered_cell_barcodes, cell_barcodes)
}

cell_covariate_matrix <- tibble(ordered_cell_barcodes = ordered_cell_barcodes, batch = batch)
write.fst(x = cell_covariate_matrix, path = paste0(processed_dir, "/cell_covariate_matrix.fst"))
saveRDS(object = genes_in_use, file = paste0(processed_dir, "/ordered_genes.RDS"))

###############
# gRNA UMI data
###############

# Load the gRNA identification information; we extract information on all guides targeting ARL15.
enh_targets_df <- read.xlsx(xlsxFile = paste0(raw_data_dir, "/enh_targets.xlsx"), sheet = 1)
arl15_region <- filter(enh_targets_df, gene_names == "ARL15") %>% pull(region)
guide_seqs <- read.xlsx(xlsxFile = paste0(raw_data_dir, "/all_oligos.xlsx"), sheet = 1) %>% rename(hg38_enh_region = "region.pos.(hg38)")
arl15_gRNA_spacer_seqs <- filter(guide_seqs, hg38_enh_region == arl15_region) %>% pull(spacer.sequence)

# Next, we create a data frame containing UMI counts for these gRNAs
raw_fs <- list.files(raw_data_dir)
gRNA_files <- paste0(raw_data_dir, "/", grep(pattern = "sgRNA-enrichment_5K-sgRNAs_Batch", x = raw_fs, value = TRUE))
res <- map(.x = gRNA_files, .f = function(curr_file) {
  print(paste("Working on file", curr_file))
  curr_gRNA_count_matrix <- read_tsv(file = curr_file, col_names = c("cell_barcode", "total_read_count", "total_umi_count", "gRNA_spacer_seqs", "read_counts", "umi_counts"), col_types = c("cccccc")) %>% select(cell_barcode, gRNA_spacer_seqs, umi_counts)
  cell_barcodes <- pull(curr_gRNA_count_matrix, cell_barcode)
  curr_batch_gRNA_umi_counts <- sapply(X = 1:nrow(curr_gRNA_count_matrix), FUN = function(row_id) {
    r <- curr_gRNA_count_matrix[row_id,]
    spacers <- str_split(r$gRNA_spacer_seqs, pattern = ";") %>% unlist()
    umi_counts <- str_split(r$umi_counts, pattern = ";") %>% unlist() %>% as.integer()
    umi_locs <- match(x = arl15_gRNA_spacer_seqs, table = spacers)
    curr_counts <- sapply(umi_locs, function(i) if (is.na(i)) 0 else umi_counts[i])
    names(curr_counts) <- arl15_gRNA_spacer_seqs
    curr_counts
  }) %>% t()
  list(cell_barcodes = cell_barcodes, umi_count_matrix = curr_batch_gRNA_umi_counts)
})

gRNA_count_matrix <- map(.x = res, .f = function(x) x$umi_count_matrix) %>% reduce(.f = rbind)
cell_barcodes <- map(.x = res, .f = function(x) x$cell_barcodes) %>% reduce(.f = c)

# Now, we modify the gRNA UMI count matrix to threshold the data; we use the method proposed by xie;
gRNA_count_matrix_thresh <- apply(X = gRNA_count_matrix, MARGIN = 2, FUN = function(column) {
  v <- sum(column >= 1)
  U <- sum(column)
  column/U > 1/v
})
gRNA_indic <- apply(X = gRNA_count_matrix_thresh, MARGIN = 1, FUN = function(r) any(r))

# Finally, confirm that the cell barcode order for the gRNA indicator matrix matches that of the cell-by-gene expression matrix and cell-specific covariate matrix.
cell_covariate_matrix <- read.fst(paste0(processed_dir, "/cell_covariate_matrix.fst"))
cell_barcodes_to_check <- pull(cell_covariate_matrix, ordered_cell_barcodes) %>% gsub(pattern = "-1", replacement = "")
m <- match(x = cell_barcodes_to_check, table = cell_barcodes)
gRNA_indic_ordered <- gRNA_indic[m]

# Put into data frame form
gRNA_indic_matrix <- data.frame(arl15_enh = gRNA_indic_ordered)
write.fst(x = gRNA_indic_matrix, path = paste0(processed_dir, "/gRNA_indicator_matrix.fst"))

##############
# Bulk RNA-seq
##############
bulk_info <- read.xlsx(xlsxFile = paste0(raw_data_dir, "/bulk_rna_info.xlsx"), sheet = 3)
bulk_info <- slice(bulk_info, 1:25) %>% select(library_name = Library.Name, gRNA = sgRNA, region = Region, biological_duplicate = Biological.Duplicate)
bulk_df <- read_tsv(file = paste0(raw_data_dir, "/GSE129825_Libraries.FeatureCounts.ARL15_enhancer.txt"), col_types = "ccccciiiiiiiiiiiiiiiiiiiiiiiiii") %>% rename("PZ788" = "PZ778", "PZ778" = "PZ778_1")
bulk_df <- filter(bulk_df, Geneid %in% all_protein_coding_genes)

write.fst(x = bulk_info, path = paste0(processed_dir, "/bulk_RNAseq_info.fst"))
write.fst(x = bulk_df, path = paste0(processed_dir, "/bulk_RNAseq.fst"))

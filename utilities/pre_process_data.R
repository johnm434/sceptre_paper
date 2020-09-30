# Pre-process the data

args <- commandArgs(trailingOnly = TRUE)
offsite_dir <- args[1] # "/Volumes/tims_new_drive/research/sceptre_files"
processed_dir <- paste0(offsite_dir, "/data/gasperini/processed")

#######################
# Gasperini CRISPR data
#######################

# 1. Create file-backed matrices for the .mtx expression matrix; both the gene-cell matrix and its transpose.
path_to_mtx <- paste0(offsite_dir, "/data/gasperini/raw/GSE120861_at_scale_screen.exprs.mtx")
dest_folder <- processed_dir
fbms <- create_fbm_from_mtx(path_to_mtx, dest_folder)
# save the .bk metadata for both of these matrices (nrow, ncell, data type, file location)
for (i in 1:length(fbms)) {
  file_name <- paste0(names(fbms)[[i]], ".rds")
  full_file_path <- paste0(dest_folder, "/",file_name)
  saveRDS(fbms[[i]], full_file_path)
}

# 2. Read in the gRNA expression data
rm(list = ls()); gc()
library(monocle, quietly = TRUE) # monocole required to load monocole CellDataSet object; we will extract only the cell-specific metadata
m_object <- readRDS("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/raw/GSE120861_at_scale_screen.cds.rds")
cell_metadata <- pData(m_object)
rm(m_object)
covariates_cols <- 1:18
gRNA_cols <- 19:ncol(cell_metadata)

gRNA_indicators <- cell_metadata[,gRNA_cols]
cell_covariates <- cell_metadata[,covariates_cols]

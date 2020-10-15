# Download data
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers_xie/paths_to_dirs.R"))
suppressPackageStartupMessages(library(R.utils))

################################
# 1. Mosaic-seq single-cell data
################################
# Set the source and the destination; perform download and untar
dest <- paste0(raw_data_dir, "/GSE129837_RAW.tar")
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129837&format=file", destfile = dest)
untar(dest, exdir = raw_data_dir)
file.remove(dest)

############################
# 2. Bulk RNA-seq validation
############################
dest <- paste0(raw_data_dir, "/GSE129825_Libraries.FeatureCounts.ARL15_enhancer.txt.gz")
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129825&format=file&file=GSE129825%5FLibraries%2EFeatureCounts%2EARL15%5Fenhancer%2Etxt%2Egz", destfile = dest)

dest <- paste0(raw_data_dir, "/GSE129825_Libraries.FeatureCounts.MYB_enhancer.txt.gz")
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129825&format=file&file=GSE129825%5FLibraries%2EFeatureCounts%2EMYB%5Fenhancer%2Etxt%2Egz", destfile = dest)

# Unzip all the txt.gz files
all_file_names <- list.files(raw_data_dir)
to_unzip <- paste0(raw_data_dir, "/", grep(pattern = '*.gz', x = all_file_names, value = TRUE))
for (file in to_unzip) {
  if (file.exists(file)) {
    gunzip(file)
  }
}

############################
# 3. Spreadsheets and tables
############################
dest <- paste0(raw_data_dir, "/all_oligos.xlsx")
download.file(url = "https://ars.els-cdn.com/content/image/1-s2.0-S2211124719313956-mmc2.xlsx", destfile = dest)

dest <- paste0(raw_data_dir, "/enh_targets.xlsx")
download.file(url = "https://ars.els-cdn.com/content/image/1-s2.0-S2211124719313956-mmc4.xlsx", destfile = dest)

dest <- paste0(raw_data_dir, "/bulk_rna_info.xlsx")
download.file(url = "https://ars.els-cdn.com/content/image/1-s2.0-S2211124719313956-mmc3.xlsx", destfile = dest)

# Download the Genes.xlsx file from this url: https://osf.io/edjzv/ (source: Human protein-coding genes and gene feature statistics in 2019 by Piovesan et al in BMC Research Notes)

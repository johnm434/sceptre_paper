####################################################################
# Download raw data from the web
# Note: The GeneHancer database is proprietary and therefore
# must be accessed piecemeal via https://genealacart.genecards.org/.
####################################################################

library(R.utils, quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers/file_paths_to_dirs.R"))

################################
# Download Gasperini et al. data
################################

# path to store the raw Gasperini data
raw_data_dir_gasp <- raw_data_dir

# URL of data
remote <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE120861&format=file&file="

# Gasperini et al results
all_deg_results_filename <- "GSE120861_all_deg_results.at_scale.txt"

# names of genes
genes_filename <- "GSE120861_at_scale_screen.genes.txt"

# names of cells
cells_filename <- "GSE120861_at_scale_screen.cells.txt"

# "reference cells" used by Gasperini et al for computational purposes
reference_cells_filename <- "GSE120861_50k_reference_cells.rds"

# all (gRNA, gene) pairs
gRNAgroup_pair_table_filename <- "GSE120861_gene_gRNAgroup_pair_table.at_scale.txt"

# list of gRNA groups used
gRNA_groups_filename <- "GSE120861_grna_groups.at_scale.txt"

# Monocle Cell Data Set object with all data
cds_filename <- "GSE120861_at_scale_screen.cds.rds"

# Expression data
expression_filename <- "GSE120861_at_scale_screen.exprs.mtx"

# list of files to download
filenames <- c(all_deg_results_filename,
              genes_filename,
              cells_filename,
              reference_cells_filename,
              cds_filename,
              expression_filename,
              gRNAgroup_pair_table_filename,
              gRNA_groups_filename)

# download files if not already present
for (filename in filenames) {
  if (!file.exists(paste0(raw_data_dir_gasp, "/", filename))) {
    cat(paste0("Downloading ", filename, "\n"))
    source <- paste0(remote, filename, ".gz")
    dest <- paste0(raw_data_dir_gasp, "/", filename, ".gz")
    download.file(source, dest)
    gunzip(paste0(dest))
  }
}

# Download supplementary Table S2 from Cell
supplementary_table_file <- "https://www.cell.com/cms/10.1016/j.cell.2018.11.029/attachment/7319ccb0-a8c0-45f3-8203-26b9159b0102/mmc2.xlsx"
download.file(supplementary_table_file, paste0(raw_data_dir_gasp, "/Gasperini_TableS2.xlsx"))

####################
# Download HI-C data
####################

# URL of data
functional_data_dir <- paste0(offsite_dir, "/data/functional")
remote <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63525&format=file&file="
TADs_filename <- "GSE63525_K562_Arrowhead_domainlist.txt"
contact_matrices_dirname <- "GSE63525_K562_intrachromosomal_contact_matrices"

# download TADs file, if necessary
dest <- paste0(functional_data_dir, "/HIC/", TADs_filename)
if (!file.exists(dest)) {
  cat(paste0("Downloading ",TADs_filename, "\n"))
  download.file(paste0(remote, TADs_filename, ".gz"),
                paste0(dest, ".gz"))
  gunzip(paste0(dest, ".gz"))
}

# download contact matrices, if necessary
dest <- paste0(functional_data_dir, "/HIC/", contact_matrices_dirname)
if (!dir.exists(dest)) {
  cat(paste0("Downloading ", contact_matrices_dirname, "\n"))
  download.file(paste0(remote, contact_matrices_dirname, ".tar.gz"),
                paste0(dest, ".tar.gz"))
  untar(paste0(dest, ".tar.gz"))
}

########################
# Download ChIP-seq data
########################

# Available courtesy of Shendure Lab at
# https://drive.google.com/drive/folders/177djZEEPV-udBkdqOtdjjV061O7s8dKj;
# original data can be downloaded at https://www.encodeproject.org/.

##########################
# Download GeneHancer data
##########################

# Available online at https://www.genecards.org/

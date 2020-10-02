args <- commandArgs(trailingOnly = TRUE)
offsite_dir <- args[1] # offsite_dir <- "/Volumes/tims_new_drive/research/sceptre_files"

# Load the gRNAs and genes to analyze
gRNA_gene_pairs <- read.fst(paste0(offsite_dir, "/data/gasperini/processed/gene_gRNA_pairs_to_study.fst"))


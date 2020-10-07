# Load sceptre
suppressPackageStartupMessages(library(sceptre))

# directory file paths
offsite_dir <- if (is.na(args[2])) "/Volumes/tims_new_drive/research/sceptre_files" else args[2]

# Define variables for fps to common directories.
processed_dir <- paste0(offsite_dir, "/data/gasperini/processed")
gene_precomp_dir <- paste0(offsite_dir, "/data/gasperini/precomp/gene")
gRNA_precomp_dir <- paste0(offsite_dir, "/data/gasperini/precomp/gRNA")
results_dir <- paste0(offsite_dir, "/results/gasperini")
raw_data_dir <- paste0(offsite_dir, "/data/gasperini/raw")
log_dir <- paste0(offsite_dir, "/logs/gasperini")

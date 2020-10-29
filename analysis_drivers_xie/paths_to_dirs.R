# Load sceptre
suppressPackageStartupMessages(library(sceptre))

# directory file paths
offsite_dir <- if (is.na(args[2])) "/Volumes/tims_new_drive/research/sceptre_files" else args[2]

# Define variables for fps to common directories.
processed_dir <- paste0(offsite_dir, "/data/xie/processed")
gene_precomp_dir <- paste0(offsite_dir, "/data/xie/precomp/gene")
gRNA_precomp_dir <- paste0(offsite_dir, "/data/xie/precomp/gRNA")
results_dir <- paste0(offsite_dir, "/results/xie/sceptre")
raw_data_dir <- paste0(offsite_dir, "/data/xie/raw")
log_dir <- paste0(offsite_dir, "/logs/xie")
results_dir_negative_binomial <-  paste0(offsite_dir, "/results/xie/negative_binomial")

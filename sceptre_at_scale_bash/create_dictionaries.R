args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
# Define the code and offsite dirs
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
param_file <- if(is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/analysis_drivers_xie/sceptre_function_args.R" else args[2]
# source the function arguments
source(param_file)

dicts <- create_and_store_dictionaries(gRNA_gene_pairs = gRNA_gene_pairs, gene_precomp_dir = gene_precomp_dir, gRNA_precomp_dir = gRNA_precomp_dir, results_dir = results_dir, pod_sizes = pod_sizes)
paste(dicts$n_pods[["gene"]], dicts$n_pods[["gRNA"]], dicts$n_pods[["pairs"]]) %>% cat

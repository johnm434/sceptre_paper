suppressPackageStartupMessages(library(sceptre))
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers/file_paths_to_dirs.R"))

# Load the gRNAs and genes to analyze
gRNA_gene_pairs <- read.fst(paste0(processed_dir, "/gene_gRNA_pairs_to_study.fst"))
dicts <- create_and_store_dictionaries(gRNA_gene_pairs = gRNA_gene_pairs, gene_precomp_dir = gene_precomp_dir, gRNA_precomp_dir = gRNA_precomp_dir, results_dir = results_dir, pod_sizes = c(gene = 200, gRNA = 1000, pair = 12000))

# Print to the standard output n_pods for genes, gRNAs, and pairs (in that order) so that the bash file knows how many pod_ids to iterate over.
paste(dicts$n_pods[["gene"]], dicts$n_pods[["gRNA"]], dicts$n_pods[["pairs"]]) %>% cat

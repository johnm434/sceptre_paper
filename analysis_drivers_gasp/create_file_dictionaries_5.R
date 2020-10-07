args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers_gasp/file_paths_to_dirs.R"))
source(paste0(code_dir, "/analysis_drivers_gasp/sceptre_function_args.R"))

# gRNA_gene_pairs <- slice(gRNA_gene_pairs, sample(1:nrow(gRNA_gene_pairs), 20, FALSE))
dicts <- create_and_store_dictionaries(gRNA_gene_pairs = gRNA_gene_pairs, gene_precomp_dir = gene_precomp_dir, gRNA_precomp_dir = gRNA_precomp_dir, results_dir = results_dir, pod_sizes = c(gene = 5, gRNA = 5, pair = 5))

# Print to the standard output n_pods for genes, gRNAs, and pairs (in that order) so that the bash file knows how many pod_ids to iterate over.
paste(dicts$n_pods[["gene"]], dicts$n_pods[["gRNA"]], dicts$n_pods[["pairs"]]) %>% cat

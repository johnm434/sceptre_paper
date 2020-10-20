args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(sceptre))
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
param_file <- if(is.na(args[2])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/analysis_drivers_xie/sceptre_function_args.R" else args[2]
source(param_file)

# Create a dictionary in the results directory
pod_id <- create_dictionary(gRNA_gene_pairs$gRNA_id, pod_sizes[["pair"]]) %>% pull(pod_id)
dictionary <- gRNA_gene_pairs %>% mutate(pod_id = pod_id)
write.fst(x = dictionary, path = paste0(offsite_dir, "/results/xie/dictionary.fst"))

# Print the number of pairs.
cat(dictionary$pod_id[nrow(dictionary)])

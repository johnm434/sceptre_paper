# aggregate results
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers/file_paths_to_dirs.R"))

x <- collect_results(results_dir)
write.fst(x, paste0(results_dir, "/all_results_size_unknown.fst"))

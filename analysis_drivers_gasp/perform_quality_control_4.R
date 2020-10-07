args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
source(paste0(code_dir, "/analysis_drivers/file_paths_to_dirs.R"))

# We perform a basic quality control for the Gasperini data; in particular, we restrict the cells in our analysis to those with at least 1 gRNA.
gRNA_indicator_matrix <- read.fst(paste0(processed_dir, "/gRNA_indicators.fst"))

# Determine which cells have no gRNA.
gRNA_indicator_counts <- apply(gRNA_indicator_matrix, 1, sum)
cells_to_keep <- which(gRNA_indicator_counts >= 1) %>% as.integer()

# Save the result
saveRDS(cells_to_keep, file = paste0(processed_dir, "/cells_to_keep.rds"))

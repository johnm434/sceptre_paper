args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
library(katsevich2020)
library(gridExtra)
source(paste0(code_dir, "/analysis_drivers_xie/paths_to_dirs.R"))


p1 <- readRDS(paste0(offsite_dir, "/figures/arl15_qqplot.rds"))
p2 <- readRDS(paste0(offsite_dir, "/figures/bulk_vs_sc_pvals.rds"))

p_all <- grid.arrange(p1, p2, nrow = 1, widths = 4:3)
ggsave(filename = paste0(offsite_dir, "/figures/combined_arl15enh_plot.pdf"), plot = p_all, width = 7.5, height = 3)

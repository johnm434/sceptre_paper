args <- commandArgs(trailingOnly = TRUE)
require(sceptre)
require(gap)
# Define the code and offsite dirs
offsite_dir <- if (is.na(args[1])) "/Volumes/tims_new_drive/research/sceptre_files" else args[1]
p_val_df <- paste0(offsite_dir, "/results/xie/all_results.fst") %>% read.fst()

pdf(file = paste0(offsite_dir, "/figures/histogram.pdf"), width = 5, height = 4)
hist(p_val_df$p_value, freq = FALSE, xlab = "p values", main = "Histogram")
abline(h = 1, col = "darkred")
dev.off()

pdf(file = paste0(offsite_dir, "/figures/qq_untransformed.pdf"), width = 5, height = 4)
qqunif(u = p_val_df$p_value, logscale = FALSE, type = "unif", main = "QQ-plot, untransformed")
dev.off()

pdf(file = paste0(offsite_dir, "/figures/qq_transformed.pdf"), width = 5, height = 4)
qqunif(u = p_val_df$p_value, logscale = TRUE, type = "unif", main = "QQ-plot, transformed")
dev.off()

p_val_df %>% mutate(p_val_adj = p.adjust(p_value, method = "BH")) %>% arrange(p_val_adj)

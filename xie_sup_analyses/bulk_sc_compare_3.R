args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
library(katsevich2020)
library(scales)
source(paste0(code_dir, "/analysis_drivers_xie/paths_to_dirs.R"))
source(paste0(code_dir, "/xie_sup_analyses/aux_objects.R"))

p_vals_sceptre <- paste0(results_dir, "/all_results.fst") %>% read.fst() %>% as_tibble()
p_vals_bulk <- paste0(offsite_dir, "/results/xie/bulk_rna_seq/pvals_arl15_enh.rds") %>% readRDS() %>% as_tibble()
p_vals_sceptre <- p_vals_sceptre %>% select(-gRNA_id) %>% mutate(p_value_adj = p.adjust(p_value, method = "BH"))

to_plot <- inner_join(p_vals_bulk, p_vals_sceptre, by = "gene_id") %>% rename(bulk_pval_adj = p_value_adj.x, bulk_pval = p_value.x, sceptre_pval_adj = p_value_adj.y, sceptre_pval = p_value.y) %>% mutate(is_arl15 = (gene_id == "ARL15"))

p <- ggplot(data = to_plot, mapping = aes(x = bulk_pval, y = sceptre_pval, col = is_arl15)) +
  geom_point(alpha = 0.5) +
  scale_colour_manual(values = c("grey60", "firebrick3")) +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  theme_bw() + theme(text = element_text(size = 12),
                     legend.background = element_rect(fill = "transparent", color = NA),
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(), legend.position = "none") +
                       xlab("Bulk RNA-seq p-value") + ylab("SCEPTRE p-value") +
  geom_point(mapping = aes(x = bulk_pval, y = sceptre_pval), data = filter(to_plot, is_arl15), size = 2) +
  annotate(geom = "text", x = 1e-11, y = 3e-7, label = "ARL15", col = "firebrick3")
saveRDS(object = p, file = paste0(offsite_dir, "/figures/bulk_vs_sc_pvals.rds"))
ggsave(filename = paste0(offsite_dir, "/figures/bulk_vs_sc_pvals.pdf"), plot = p, scale = 1, width = 3.5, height = 3)

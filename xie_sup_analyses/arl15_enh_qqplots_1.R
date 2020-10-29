args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE/sceptre_paper/" else args[1]
require(katsevich2020)
require(scales)
source(paste0(code_dir, "/analysis_drivers_xie/paths_to_dirs.R"))
source(paste0(code_dir, "/xie_sup_analyses/aux_objects.R"))

p_vals_sceptre <- paste0(results_dir, "/all_results.fst") %>% read.fst()
p_vals_nb <- paste0(results_dir_negative_binomial) %>% list.files(full.names = TRUE) %>% map(read.fst) %>% reduce(rbind)
p_vals_hypergeo <- paste0(processed_dir, "/hypergeometric_arl15enh_pvals_down.rds") %>% readRDS()
p_vals_hypergeo <- p_vals_hypergeo[names(p_vals_hypergeo) %in% p_vals_sceptre$gene_id]
to_plot <- tibble(method = rep(x = c("SCEPTRE", "Improved NB", "Hypergeometric"), each = length(p_vals_hypergeo)) %>% factor(), p_value = c(p_vals_sceptre$p_value, p_vals_nb$p_value, p_vals_hypergeo), gene = c(p_vals_sceptre$gene_id, p_vals_nb$gene_id, names(p_vals_hypergeo)))

ci <- 0.95
truncate_thresh <- 1e-9
qq_data <- to_plot %>%
  rename(pvalue = p_value) %>%
  group_by(method) %>%
  mutate(r = rank(pvalue), expected = ppoints(n())[r],
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
  ungroup() %>%
  mutate(pvalue = ifelse(pvalue < truncate_thresh, truncate_thresh, pvalue))

annotation_df <- filter(qq_data, gene == "ARL15")
arrow_coords <- tibble(x1 = 4e-4, x2 = annotation_df$expected + 1e-5, y1 = 1e-8, y2 = annotation_df$pvalue)

p <- qq_data %>%
  ggplot(aes(x = expected, y = pvalue, group = method, ymin = clower, ymax = cupper)) +
  geom_point(aes(color = method), size = 1, alpha = 0.5) +
  geom_ribbon(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_manual(values = setNames(plot_colors[c("hypergeometric", "improved_nb", "sceptre")], NULL), name = "Method") +
  xlab("Expected null p-value") +
  ylab("Observed p-value") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  theme_bw() + theme(legend.position = c(0.25,0.8), text = element_text(size = 12),
                     legend.background = element_rect(fill = "transparent", colour = NA),
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line()) +
                     annotate(geom = "text", x = 10e-4, y = 1e-8, label = "ARL15", col = "firebrick3") +
                     geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), lwd = 0.2, data = arrow_coords, arrow = arrow(length=unit(0.2,"cm")), inherit.aes = FALSE, col = "grey40")
saveRDS(object = p, file = paste0(offsite_dir, "/figures/arl15_qqplot.rds"))
ggsave(filename = paste0(offsite_dir, "/figures/arl15_enh_qqplot.pdf"), plot = p, scale = 1, width = 4.5, height = 3)

pdf(file = paste0(offsite_dir, "/figures/sceptre_vs_nb_pvals_arl15enh.pdf"), width = 4, height = 3)
plot(x = -log(p_vals_sceptre$p_value),y =  -log(p_vals_nb$p_value), xlab = "SCEPTRE p-values", ylab = "NB p-values")
dev.off()

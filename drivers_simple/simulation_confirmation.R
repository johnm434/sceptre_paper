require(sceptre)

n_rep <- 100
res <- sapply(1:n_rep, function(i) {
  cat(paste0("Replicate ", i, "\n"))
  sim_dat <- simulate_crispr_screen_data(num_cells = 1000,
                                         grna_mean_prob = 0.2,
                                         covariate_sampler = list(cell_size = rnorm, cell_cycle = runif),
                                         mRNA_mean_expression = 40,
                                         gRNA_effect = 0,
                                         covariate_effects = c(0.5, 1),
                                         zero_inflation = 0,
                                         neg_binom_size = 2)
  expressions <- sim_dat$Y
  gRNA_indicators <- sim_dat$X
  covariate_matrix <-sim_dat$covariate_df
  run_sceptre_gRNA_gene_pair(expressions = expressions, gRNA_indicators = gRNA_indicators, covariate_matrix = covariate_matrix)
})
# simulate data
qqplot(res, runif(1000), pch = 19, cex = 0.5)
abline(a = 0, b = 1, col = "red")
ks.test(res, punif)
hist(res, main = "p-values under null", xlab = "")
saveRDS(res, "null_p_vals.rds")



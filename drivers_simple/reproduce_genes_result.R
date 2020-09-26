# reproduce result
# The goal of this analysis is to reproduce the result obtained by Gene on a random gRNA-gene pair.
require(sceptre)
load("/Users/timbarry/Box/SCEPTRE/sceptre_paper/drivers_simple/single_pair.Rda")

# set the expression vector, gRNA indicator vector, and covariate matrix.
expressions <- df$gene_exp
gRNA_indicators <- df$grna_group_indicator
covariate_matrix <- summarize(df, percent.mito = percent.mito, prep_batch = factor(prep_batch), total_umis = log(total_umis), guide_count = log(guide_count), gene_count = log(gene_count))
result <- run_sceptre_gRNA_gene_pair(expressions = expressions, gRNA_indicators = gRNA_indicators, covariate_matrix = covariate_matrix, gene_precomp_dispersion = size, seed = 1234)

# check if the results match
abs(result$p_value_raw - random_pair_results$corrected_pvalue_raw) < 1e-5
abs(result$p_value_skew_t - random_pair_results$corrected_pvalue_st) < 1e-5
abs(result$z_value_star - random_pair_results$original_zvalue) < 1e-5

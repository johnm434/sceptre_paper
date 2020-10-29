#' Simulate CRISPR screen data
#'
#' @param n_cells number of cells
#' @param grna_mean_prob mean probability of gRNA presence across enhancers (absent covariates)
#' @param mRNA_mean_expression mean expression level of gene across genes (absent covariates)
#' @param covariate_effects_gRNA the effect of the technical covariate on gRNA presence probability; length equals number of enhancers
#' @param covariate_effects_gene the effect of the technical covariate on mean gene expression; length equals number of genes
#' @param zero_inflation zero inflation rate (between 0 and 1)
#' @param neg_binom_size size of the negative binomial distributions for gene expression
#' @param seed (optonal) seed to the random number generator
#'
#' @return list with components (i) cell_by_gene_expressions, (ii) cell_by_enhancer_perturbation_indicators, and (iii) Z
#' @export
#'
#' @examples
#' r <- simulate_crispr_screen_data(n_cells = 1000,
#' grna_mean_prob = 0.2,
#' mRNA_mean_expression = 40,
#' covariate_effects_gRNA = rep(4, 10),
#' covariate_effects_gene = rep(2, 50),
#' zero_inflation = 0,
#' neg_binom_size = 2
#' )

simulate_crispr_screen_data <- function(n_cells, grna_mean_prob, mRNA_mean_expression, covariate_effects_gRNA, covariate_effects_gene, zero_inflation, neg_binom_size, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Z <- rnorm(n_cells)
  cell_by_enhancer_perturbation_prob <- sapply(X = covariate_effects_gRNA, FUN = function(covariate_effect) {
    1 / (1 + exp(-(binomial()$linkfun(grna_mean_prob) + covariate_effect * Z)))
  })
  cell_by_enhancer_perturbation_indicators <- apply(X = cell_by_enhancer_perturbation_prob, MARGIN = 2, FUN = function(col) {
    rbinom(n = n_cells, size = 1, prob = col)
  })
  cell_by_gene_mean_expression <- sapply(X = covariate_effects_gene, FUN = function(covariate_effect) {
    exp(log(mRNA_mean_expression) + covariate_effect * Z)
  })
  cell_by_gene_expressions <- apply(cell_by_gene_mean_expression, MARGIN = 2, FUN = function(col) {
    rnbinom(n = n_cells, size = neg_binom_size, mu = col)
  })
  n_genes <- length(covariate_effects_gene)
  cell_by_gene_expressions <- cell_by_gene_expressions * matrix(data = rbinom(n = n_genes * n_cells, size = 1, prob = 1 - zero_inflation), nrow = n_cells, ncol = n_genes)
  return(list(cell_by_gene_expressions = cell_by_gene_expressions, cell_by_enhancer_perturbation_indicators = cell_by_enhancer_perturbation_indicators, Z = Z))
}

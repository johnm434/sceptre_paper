#' Simulate CRISPR screen data.
#'
#' A function to simulate pooled CRISPR screen data.
#'
#' @param num_cells number of cells in experiment
#' @param grna_mean_prob mean number of cells perturbed (absent covariates)
#' @param covariate_sampler a list of functions. Each function in this list takes as an argument num_cells and returns num_cells random variates.
#' @param mRNA_mean_expression mean mRNA expression of the cell (absent covariates)
#' @param gRNA_effect effect of the gRNA perturbation on mRNA expression
#' @param covariate_effects a numeric vector indicating the covariate effects on mRNA expression
#' @param zero_inflation a numeric scalar between 0 and 1 (inclusive) giving the mean fraction of cell expressions set to zero
#' @param neg_binom_size size parameter for the negative binomial model
#'
#' @return
#' a list containing (i) the expression vector Y, (ii) the data frame of covariates, and (iii) the gRNA presence indicator vector.
#' @export
#'
#' @examples
#' num_cells <- 1000
#' grna_mean_prob <- 0.2
#' covariate_sampler <- list(
#'  cell_size = runif,
#'  cell_cycle = function(x) {runif(n = x, min = 0, max = 1)}
#' )
#' mRNA_mean_expression <- 40
#' gRNA_effect <- 4
#' covariate_effects <- c(2, 1)
#' zero_inflation <- 0
#' neg_binom_size <- 2
#' simulated_data <- simulate_crispr_screen_data(num_cells, grna_mean_prob, covariate_sampler, mRNA_mean_expression, gRNA_effect, covariate_effects, zero_inflation, neg_binom_size)

simulate_crispr_screen_data <- function(num_cells, grna_mean_prob, covariate_sampler, mRNA_mean_expression, gRNA_effect, covariate_effects, zero_inflation, neg_binom_size) {
  dat <- map_dfr(.x = covariate_sampler, function(f) f(num_cells))
  formula <- paste0("~", paste0(colnames(dat), collapse = " + ")) %>% as.formula
  cov_model <- model.matrix(formula, data = dat)
  l_i_g <- as.numeric(cov_model %*% c(logitlink(grna_mean_prob), covariate_effects))
  pi_true <-  1/(1 + exp(-(l_i_g)))
  X <- rbinom(n = num_cells, size = 1, prob = pi_true)
  updated_cov_model <- cbind(cov_model, X)
  l_i_m <- as.numeric(updated_cov_model %*% c(log(mRNA_mean_expression), covariate_effects, gRNA_effect))
  Y <- rnbinom(n = num_cells,
              size = neg_binom_size,
              mu = exp(l_i_m))
  Y <- Y * rbinom(n = num_cells,
                 size = 1,
                 prob = 1 - zero_inflation)
  return(list(Y = Y, covariate_matrix = dat, X = X))
}

#' Run sceptre using precomputations for gRNAs and genes.
#'
#' This function is the workhorse function of the sceptre package. It runs a distilled CRT using a negative binomial test statistic based on an expression vector, a gRNA indicator vector, an offset vector (from the distillation step), gRNA conditional probabilities, an estimate of the negative binomial dispersion parameter, and the number of resampling replicates.
#'
#' This is a one-tailed, left-sided test. Activators.
#' @param expressions a vector of gene expressions (in UMI counts)
#' @param gRNA_indicators a vector of gRNA indicators
#' @param gRNA_precomp a vector of conditional probabilities for gRNA assignments
#' @param gene_precomp_dispersion the pre-computed dispersion
#' @param gene_precomp_offsets the pre-computed distillation offsets
#' @param B the number of resamples to make (default 500)
#' @param seed an arguement to set.seed; if null, no seed is set
#'
#' @export
#' @return a p-value of the null hypothesis of no gRNA effect on gene expression
run_sceptre_using_precomp <- function(expressions, gRNA_indicators, gRNA_precomp, gene_precomp_dispersion, gene_precomp_offsets, B, seed) {
  if (!is.null(seed)) set.seed(seed)

  # compute the test statistic on the real data
  fit_star <- vglm(formula = expressions[gRNA_indicators == 1] ~ 1, family = negbinomial.size(gene_precomp_dispersion), offset = gene_precomp_offsets[gRNA_indicators == 1])
  t_star <- summaryvglm(fit_star)@coef3["(Intercept)", "z value"]

  # resample B times
  t_nulls <- sapply(1:B, function(i) {
    if (i %% 100 == 0) cat(paste0("Running resample ", i ,"/", B, "\n"))
    gRNA_indicators_null <- rbinom(n = length(gRNA_precomp), size = 1, prob = gRNA_precomp)
    fit_null <- vglm(formula = expressions[gRNA_indicators_null == 1] ~ 1, family = negbinomial.size(gene_precomp_dispersion), offset = gene_precomp_offsets[gRNA_indicators_null == 1])
    summaryvglm(fit_null)@coef3["(Intercept)", "z value"]
  })

  # Fit a skew-t distribution and obtain a p-value
  p_value_skew_t <- NA
  skew_t_fit <- tryCatch(selm(t_nulls ~ 1, family = "ST"), error = function(e) return(NA), warning = function(w) return(NA))
  if (class(skew_t_fit) == "selm") { # If the fit worked,
    dp <- skew_t_fit@param$dp # then extract the parameters.
    if (!any(is.na(dp))) { # If all the fitted parameters are numbers,
      p_value_skew_t <- pst(x = t_star, dp = dp) # then compute the skew t-based p-value.
    }
  }
  p_value_raw <- mean(c(-Inf, t_nulls) <= t_star)
  return(list(p_value_raw = p_value_raw, p_value_skew_t = p_value_skew_t, z_value_star = t_star))
}


#' Run sceptre on a gRNA-gene pair
#'
#' This function runs the sceptre algorithm on a gRNA-gene pair. It requires as arguments the gene expression vector, the gRNA indicator vector, and the covariate matrix. Users optionally can pass the gRNA precomputation or gene precomputation as arguments.
#'
#' If the user wants to pass both the gRNA precomputation and gene precomputation, the user should use instead the function run_sceptre_using_precomp. Most users of this pacakge will NOT pass either gRNA_precomp or gene_precomp as arguements to this function.
#'
#' @param expressions a vector a gene expressions
#' @param gRNA_indicators a vector of gRNA inicators
#' @param covariate_matrix the matrix of cell-specific covariates (e.g., library size, batch effect, cell cycle, etc.)
#' @param gRNA_precomp (optional) the gRNA precomputation (a vector of gRNA presence conditional probabilities)
#' @param gene_precomp_dispersion (optional) the pre-computed gene dispersion
#' @param gene_precomp_offsets (optional) the pre-computed gene offsets
#' @param B  number of resamples (default 500)
#'
#' @return a p-value of the null hypothesis of no gRNA effect on gene expression
#' @export
#'
#' @examples
#' # An example in which the alternative is true.
#' sim_dat <- simulate_crispr_screen_data(num_cells = 1000,
#' grna_mean_prob = 0.2,
#' covariate_sampler = list(cell_size = rnorm, cell_cycle = runif),
#' mRNA_mean_expression = 40,
#' gRNA_effect = -4,
#' covariate_effects = c(0.5, 1),
#' zero_inflation = 0,
#' neg_binom_size = 2)
#' expressions <- sim_dat$Y
#' gRNA_indicators <- sim_dat$X
#' covariate_matrix <-sim_dat$covariate_df
#' run_sceptre_gRNA_gene_pair(expressions = expressions, gRNA_indicators = gRNA_indicators, covariate_matrix = covariate_matrix)
#'
#' # An example in which the null is true.
#' sim_dat <- simulate_crispr_screen_data(num_cells = 1000,
#' grna_mean_prob = 0.2,
#' covariate_sampler = list(cell_size = rnorm, cell_cycle = runif),
#' mRNA_mean_expression = 40,
#' gRNA_effect = 0,
#' covariate_effects = c(0.5, 1),
#' zero_inflation = 0,
#' neg_binom_size = 2)
#' expressions <- sim_dat$Y
#' gRNA_indicators <- sim_dat$X
#' covariate_matrix <-sim_dat$covariate_df
#' run_sceptre_gRNA_gene_pair(expressions = expressions, gRNA_indicators = gRNA_indicators, covariate_matrix = covariate_matrix)
run_sceptre_gRNA_gene_pair <- function(expressions, gRNA_indicators, covariate_matrix, gRNA_precomp = NULL, gene_precomp_dispersion = NULL, gene_precomp_offsets = NULL, B = 500, seed = NULL) {
  if (is.null(gRNA_precomp)) gRNA_precomp <- run_gRNA_precomputation(gRNA_indicators, covariate_matrix)
  if (is.null(gene_precomp_dispersion) || is.null(gene_precomp_offsets)) gene_precomp <- run_gene_precomputation(expressions, covariate_matrix, gene_precomp_dispersion, gene_precomp_offsets)
  out <- run_sceptre_using_precomp(expressions = expressions, gRNA_indicators = gRNA_indicators, gRNA_precomp = gRNA_precomp, gene_precomp_dispersion = gene_precomp$gene_precomp_dispersion, gene_precomp_offsets = gene_precomp$gene_precomp_offsets, B = B, seed = seed)
  return(out)
}


#' Run gRNA precomputation
#'
#' This function runs the precomputation for a given gRNA.
#'
#' @param gRNA_indicators a vector of gRNA indicators
#' @param covariate_matrix the cell-specific covariate matrix
#'
#' @export
#' @return the fitted probabilities
#'
#' @examples
#' sim_dat <- simulate_crispr_screen_data(num_cells = 1000,
#' grna_mean_prob = 0.2,
#' covariate_sampler = list(cell_size = rnorm, cell_cycle = runif),
#' mRNA_mean_expression = 40,
#' gRNA_effect = -4,
#' covariate_effects = c(0.5, 1),
#' zero_inflation = 0,
#' neg_binom_size = 2)
#' gRNA_indicators <- sim_dat$X
#' covariate_matrix <-sim_dat$covariate_df
#' fitted_probs <- run_gRNA_precomputation(gRNA_indicators, covariate_matrix)
run_gRNA_precomputation <- function(gRNA_indicators, covariate_matrix) {
  fit_model_grna <- glm(gRNA_indicators ~ ., family = binomial(), data = covariate_matrix)
  out <- as.numeric(fitted(fit_model_grna))
  return(out)
}

#' Run gene precomputation
#'
#' This function runs the precomputation for a given gene. In particlar, it fits an NB regression of expression against covariates. The estimate of theta (i.e., the NB dispersion parameter) is obtained from glm.nb function. This is sensible as, under the null hypothesis, the NB model without the gRNA indicator is true. Offsets are obtained by log-transforming the fitted values.
#'
#' @param expressions the vector of gene expressions
#' @param covariate_matrix the cell-specific covariate matrix
#'
#' @return a named list containing two items: offsets and dispersion.
#' @export
#'
#' @examples
#' sim_dat <- simulate_crispr_screen_data(num_cells = 1000,
#' grna_mean_prob = 0.2,
#' covariate_sampler = list(cell_size = rnorm, cell_cycle = runif),
#' mRNA_mean_expression = 40,
#' gRNA_effect = 0,
#' covariate_effects = c(0.5, 1),
#' zero_inflation = 0,
#' neg_binom_size = 2)
#' expressions <- sim_dat$Y
#' covariate_matrix <- sim_dat$covariate_df
#' gene_precomp <- run_gene_precomputation(expressions, covariate_matrix, NULL, NULL)
run_gene_precomputation <- function(expressions, covariate_matrix, gene_precomp_dispersion, gene_precomp_offsets) {
  # cases on gene_precomp_dispersion
  if (is.null(gene_precomp_dispersion)) { # no dispersion supplied; use glm.nb to estimate dispersion and fit model
    tryCatch({
      fit_nb <- glm.nb(formula = expressions ~ ., data = covariate_matrix)
      fitted_vals <- as.numeric(fit_nb$fitted.values)
      gene_precomp_dispersion_out <- fit_nb$theta
    }, error = function(e) {
      pois_fit <- glm(expressions ~ ., data = covariate_matrix, family = poisson())
      gene_precomp_dispersion_out <- theta.ml(expressions, pois_fit$fitted.values)[1]
      fit_nb <- vglm(formula = expressions ~ ., family = negbinomial.size(gene_precomp_dispersion_out), data = covariate_matrix)
      fitted_vals <- as.numeric(fittedvlm(fit_nb))
    })
  } else { # dispersion supplied; use vglm to fit model
    gene_precomp_dispersion_out <- gene_precomp_dispersion
    fit_nb <- vglm(formula = expressions ~ ., family = negbinomial.size(gene_precomp_dispersion_out), data = covariate_matrix)
    fitted_vals <- as.numeric(fittedvlm(fit_nb))
  }
  if (is.null(gene_precomp_offsets)) {
    gene_precomp_offsets_out <- log(fitted_vals)
  } else {
    gene_precomp_offsets_out <- gene_precomp_offsets
  }
  out <- list(gene_precomp_offsets = gene_precomp_offsets_out, gene_precomp_dispersion = gene_precomp_dispersion_out)
  return(out)
}

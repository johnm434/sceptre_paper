#' Run sceptre using precomputations for gRNAs and genes.
#'
#' This function is the workhorse function of the sceptre package. It runs a distilled CRT using a negative binomial test statistic based on an expression vector, a gRNA indicator vector, an offset vector (from the distillation step), gRNA conditional probabilities, an estimate of the negative binomial size parameter, and the number of resampling replicates.
#'
#' This currently is a one-tailed, left-sided test. Thus, it is suitable for up-regulatory elements like enhancers and promoters but not down-regulatory elements like silencers.
#'
#' @param expressions a vector of gene expressions (in UMI counts)
#' @param gRNA_indicators a vector of gRNA indicators
#' @param gRNA_precomp a vector of conditional probabilities for gRNA assignments
#' @param gene_precomp_size the pre-computed size
#' @param gene_precomp_offsets the pre-computed distillation offsets
#' @param B the number of resamples to make (default 500)
#' @param seed an arguement to set.seed; if null, no seed is set
#'
#' @export
#' @return a p-value of the null hypothesis of no gRNA effect on gene expression
run_sceptre_using_precomp <- function(expressions, gRNA_indicators, gRNA_precomp, gene_precomp_size, gene_precomp_offsets, B, seed) {
  if (!is.null(seed)) set.seed(seed)

  # compute the test statistic on the real data
  fit_star <- vglm(formula = expressions[gRNA_indicators == 1] ~ 1, family = negbinomial.size(gene_precomp_size), offset = gene_precomp_offsets[gRNA_indicators == 1])
  t_star <- summaryvglm(fit_star)@coef3["(Intercept)", "z value"]

  # resample B times
  t_nulls <- sapply(1:B, function(i) {
    if (i %% 100 == 0) cat(paste0("Running resample ", i ,"/", B, ".\n"))
    gRNA_indicators_null <- rbinom(n = length(gRNA_precomp), size = 1, prob = gRNA_precomp)
    tryCatch({
      fit_null <- vglm(formula = expressions[gRNA_indicators_null == 1] ~ 1, family = negbinomial.size(gene_precomp_size), offset = gene_precomp_offsets[gRNA_indicators_null == 1])
      summaryvglm(fit_null)@coef3["(Intercept)", "z value"]},
      error = function(e) return(NA),
      warning = function(w) return(NA)
      )
  })
  t_nulls <- t_nulls[!is.na(t_nulls)]

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
  out <- if (is.na(p_value_skew_t)) p_value_raw else p_value_skew_t
  return(out)
}


#' Run sceptre on a gRNA-gene pair
#'
#' This function runs the sceptre algorithm on a single gRNA-gene pair. It requires as arguments the gene expression vector, the gRNA indicator vector, and the covariate matrix. Users optionally can pass the gRNA precomputation or gene precomputation as arguments.
#'
#' @param expressions a vector a gene expressions
#' @param gRNA_indicators a vector of gRNA inicators
#' @param covariate_matrix the matrix of cell-specific covariates (e.g., library size, batch effect, cell cycle, etc.)
#' @param gene_precomp_size (optional) the pre-computed size of the gene NB distribution
#' @param B number of resamples (default 500)
#' @param seed (optional) seed to the random number generator
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
#' run_sceptre_gRNA_gene_pair(expressions = expressions,
#' gRNA_indicators = gRNA_indicators,
#' covariate_matrix = covariate_matrix)
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
#' run_sceptre_gRNA_gene_pair(expressions = expressions,
#' gRNA_indicators = gRNA_indicators,
#' covariate_matrix = covariate_matrix)
run_sceptre_gRNA_gene_pair <- function(expressions, gRNA_indicators, covariate_matrix, gene_precomp_size = NULL, B = 500, seed = NULL) {
  cat(paste0("Running gRNA precomputation.\n"))
  gRNA_precomp <- run_gRNA_precomputation(gRNA_indicators, covariate_matrix)

  cat(paste0("Running gene precomputation.\n"))
  gene_precomp <- run_gene_precomputation(expressions, covariate_matrix, gene_precomp_size)

  out <- run_sceptre_using_precomp(expressions = expressions, gRNA_indicators = gRNA_indicators, gRNA_precomp = gRNA_precomp, gene_precomp_size = gene_precomp$gene_precomp_size, gene_precomp_offsets = gene_precomp$gene_precomp_offsets, B = B, seed = seed)
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
#' This function runs the precomputation for a given gene. In particlar, it fits an NB regression of expression against covariates. The estimate of theta (i.e., the NB size parameter) is obtained from glm.nb function. This is sensible as, under the null hypothesis, the NB model without the gRNA indicator is true. Offsets are obtained by log-transforming the fitted values.
#'
#' @param expressions the vector of gene expressions
#' @param covariate_matrix the cell-specific covariate matrix
#' @param gene_precomp_size the pre-computed size parameter (NULL if none)
#'
#' @return a named list containing two items: offsets and size.
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
run_gene_precomputation <- function(expressions, covariate_matrix, gene_precomp_size) {
  # cases on gene_precomp_size
  if (is.null(gene_precomp_size)) { # no size supplied; use glm.nb to estimate size and fit model

    backup_2 <- function(pois_fit) {
      theta.mm(y = expressions, mu = pois_fit$fitted.values, dfr = pois_fit$df.residual)
    }
    backup <- function() {
      pois_fit <- glm(expressions ~ ., data = covariate_matrix, family = poisson())
      gene_precomp_size_out <- tryCatch({
        theta.ml(expressions, pois_fit$fitted.values, limit = 50)[1]
      }, error = function(e) backup_2(pois_fit), warning = function(w) backup_2(pois_fit))
      fit_nb <- vglm(formula = expressions ~ ., family = negbinomial.size(gene_precomp_size_out), data = covariate_matrix)
      fitted_vals <- as.numeric(fittedvlm(fit_nb))
      list(fitted_vals = fitted_vals, gene_precomp_size_out = gene_precomp_size_out)
    }

    result <- tryCatch({
      fit_nb <- glm.nb(formula = expressions ~ ., data = covariate_matrix)
      fitted_vals <- as.numeric(fit_nb$fitted.values)
      gene_precomp_size_out <- fit_nb$theta
      list(fitted_vals = fitted_vals, gene_precomp_size_out = gene_precomp_size_out)
    }, error = function(e) backup(), warning = function(w) backup())

    fitted_vals <- result$fitted_vals; gene_precomp_size_out <- result$gene_precomp_size_out

  } else { # size supplied; use vglm to fit model
    gene_precomp_size_out <- gene_precomp_size
    fit_nb <- vglm(formula = expressions ~ ., family = negbinomial.size(gene_precomp_size_out), data = covariate_matrix)
    fitted_vals <- as.numeric(fittedvlm(fit_nb))
  }

  gene_precomp_offsets_out <- log(fitted_vals)
  out <- list(gene_precomp_offsets = gene_precomp_offsets_out, gene_precomp_size = gene_precomp_size_out)
  return(out)
}

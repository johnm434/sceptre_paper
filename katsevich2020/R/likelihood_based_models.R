# This file contains functions for performing likelihood-based inference on single-cell pooled crispr screen data.


#' Run negative binomial model
#'
#' Runs a negative binomial model on a vector of expression data, gRNA indicators, and a covariate matrix.
#'
#' @param expressions a vector of expressions
#' @param gRNA_indicators a vector of gRNA indicators
#' @param covariate_matrix a covariate matrix
#'
#' @return p-value corresponding to the gRNA indicator
#' @export
#'
#' @examples
#' offsite_dir <- "/Volumes/tims_new_drive/research/sceptre_files"
#' processed_dir <- paste0(offsite_dir, "/data/xie/processed")
#' covariate_matrix <- read.fst(paste0(processed_dir, "/covariate_model_matrix.fst"))
#' ordered_gene_ids <- readRDS(paste0(processed_dir, "/ordered_genes.RDS"))
#' cell_gene_expression_matrix <- readRDS(paste0(processed_dir, "/exp_mat_metadata.rds")) %>% load_fbm
#' expressions <- cell_gene_expression_matrix[, which(ordered_gene_ids == "ARL15")]
#' gRNA_indicators <- paste0(processed_dir, "/gRNA_indicator_matrix.fst") %>% read.fst() %>% pull()
run_NB_model <- function(expressions, gRNA_indicators, covariate_matrix) {
  # construct the full model matrix
  if("gRNA_indic" %in% colnames(covariate_matrix)) stop("gRNA_indic should be passed as a vector and not included as a column in the covariate matrix.")
  model_matrix <- mutate(covariate_matrix, gRNA_indic = as.integer(gRNA_indicators))

  # first, we estimate the size parameter.
  backup_2 <- function(pois_fit) {
    theta.mm(y = expressions, mu = pois_fit$fitted.values, dfr = pois_fit$df.residual)
  }
  backup <- function() {
    pois_fit <- glm(expressions ~ ., data = model_matrix, family = poisson())
    tryCatch({
      theta.ml(y = expressions, mu = pois_fit$fitted.values, limit = 50)[1]
    }, error = function(e) backup_2(pois_fit), warning = function(w) backup_2(pois_fit))
  }
  size <- tryCatch({
  fit_nb <- glm.nb(formula = expressions ~ ., data = model_matrix)
  fit_nb$theta
  }, error = function(e) backup(), warning = function(w) backup())

  # next, we fit the model with known size parameter.
  fit_star <- vglm(formula = expressions ~ ., family = negbinomial.size(size), data = model_matrix)
  # we extract the z-value corresponding to gRNA_indic
  z_val <- summaryvglm(fit_star)@coef3["gRNA_indic", "z value"]
  # Return the appropriate one-sided p-value
  p <- pnorm(q = z_val, lower.tail = TRUE)
  return(p)
}

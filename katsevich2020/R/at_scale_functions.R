# At scale functions

#' Map reduce at scale
#'
#' @param f function we are applying
#' @param pod_id pod id
#' @param dictionary dictionary of gRNA-gene pairs
#' @param results_dir directory in which to store the result
#' @param cell_gene_expression_matrix the cell-by-gene expression matrix
#' @param ordered_gene_ids the ordered column nmes of the expression matrix
#' @param gRNA_indicator_matrix_fp a file-path to the gRNA indicator matrix
#' @param covariate_matrix matrix of covariates
#' @param cell_subset integer vector indicating cell subset
#' @param log_dir location to sink log fles
#' @export
map_reduce_at_scale <- function(f, pod_id, dictionary, results_dir, cell_gene_expression_matrix, ordered_gene_ids, gRNA_indicator_matrix_fp, covariate_matrix, cell_subset = NULL, log_dir = NULL) {
  if (!is.null(log_dir)) activate_sink(paste0(log_dir, "/" , f, "_", pod_id, ".Rout"))
  pairs_to_analyze <- filter(dictionary, pod_id == !!pod_id)
  p_vals <- sapply(X = 1:nrow(pairs_to_analyze), FUN = function(i) {
    gene_id <- pairs_to_analyze[[i, "gene_id"]]
    gRNA_id <- pairs_to_analyze[[i, "gRNA_id"]]
    cat(paste("Running", f, "on gene", gene_id, "and gRNA", gRNA_id, "\n"))
    expressions <- cell_gene_expression_matrix[, which(ordered_gene_ids == gene_id)]
    gRNA_indicators <- read.fst(gRNA_indicator_matrix_fp, columns = gRNA_id) %>% pull()
    # subset if necessary
    if (!is.null(cell_subset)) {
      expressions <- expressions[cell_subset]
      gRNA_indicators <- gRNA_indicators[cell_subset]
      covariate_matrix <- covariate_matrix[cell_subset,]
    }
    do.call(what = f, args = list(expressions = expressions, gRNA_indicators = gRNA_indicators, covariate_matrix = covariate_matrix))
  })
  out <- pairs_to_analyze %>% summarize(gene_id = gene_id, gRNA_id = gRNA_id) %>% mutate(p_value = p_vals)
  write.fst(x = out, path = paste0(results_dir, "/result_", pod_id, ".fst"))
  if (!is.null(log_dir)) deactivate_sink()
}

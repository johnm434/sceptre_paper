# Depricated

#' Run gene precomputation at scale (depricated)
#'
#' @param pod_id ID of the gene pod on which to run the precomputation
#' @param gene_precomp_dir file path to the gene precomputation directory
#' @param cell_gene_expression_matrix an FBM containing the expression data (rows cells, columns genes)
#' @param ordered_gene_ids the gene IDs (i.e., names)
#' @param covariate_matrix the cell-specific covariate matrix
#' @param cell_subset (optional) integer vector identifying the cells to use in the model
#' @param gene_sizes (optional) a vector of already-estimated (or known) gene sizes
#' @param log_dir file path to the directory in which to sink the log output
#'
#' @return NULL
#' @export
#'
#' @examples
#' pod_id <- 1
#' gene_precomp_dir <- "/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/precomp/gene/"
#' cell_gene_expression_matrix_metadata <- readRDS("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/expression_FBM_metadata.rds")
#' cell_gene_expression_matrix <- load_fbm(cell_gene_expression_matrix_metadata)
#' covariate_matrix <- read.fst("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/cell_covariate_model_matrix.fst")
#' cell_subset <- readRDS("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/cells_to_keep.rds")
#' ordered_gene_ids <- readRDS("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/ordered_gene_ids.RDS")
run_gene_precomputation_at_scale <- function(pod_id, gene_precomp_dir, cell_gene_expression_matrix, ordered_gene_ids, covariate_matrix, cell_subset = NULL, log_dir = NULL, gene_sizes = NULL) {
  if (!is.null(log_dir)) activate_sink(paste0(log_dir, "/gene_precomp_", pod_id, ".Rout"))
  # subset covariate matrix by rows
  if (!is.null(cell_subset)) covariate_matrix <- covariate_matrix[cell_subset,]
  # obtain the genes on which to preform the precomputation
  gene_dictionary <- read.fst(paste0(gene_precomp_dir, "/gene_dictionary.fst")) %>% filter(pod_id == !!pod_id)
  gene_ids <- gene_dictionary %>% pull(id) %>% as.character()

  # Run the precomputations
  precomps <- map(gene_ids, function(gene_id) {
    cat(paste0("Running precomputation for gene ", gene_id, ".\n"))
    integer_id <- which(gene_id == ordered_gene_ids)
    expressions <- cell_gene_expression_matrix[,integer_id]
    if (!is.null(cell_subset)) expressions <- expressions[cell_subset]
    curr_gene_size <- gene_sizes[[gene_id]]
    run_gene_precomputation(expressions = expressions, covariate_matrix = covariate_matrix, gene_precomp_size = curr_gene_size)
  })
  names(precomps) <- gene_ids
  out_offsets <- map_dfc(precomps, function(l) l$gene_precomp_offsets)
  out_sizes <- map_dbl(precomps, function(l) l$gene_precomp_size)

  # save the precomputations
  offset_save_fp <- (gene_dictionary %>% pull(offset_file))[1] %>% as.character()
  size_save_fp <- (gene_dictionary %>% pull(size_file))[1] %>% as.character()
  write.fst(out_offsets, offset_save_fp)
  saveRDS(out_sizes, size_save_fp)
  if (!is.null(log_dir)) deactivate_sink()
}

# scaling up package
# This script provides functions to scale up sceptre to large datasets.

#' Activate and deactivate sink
#'
#' A convenience function for programs that run at scale. activate_sink activates the sink to the log file, and deactivate_sink deactivates the sink.
#'
#' @param log_file_name the name of the log file in which to sink the output (including the file path)
#'
#' @return
activate_sink <- function(log_file_name) {
  sink(log_file_name)
  sink(stdout(), type = "message")
}

#' @rdname activate_sink
deactivate_sink <- function() {
  sink(NULL, type = "message")
  sink()
}

#' Create FBM from .mtx file
#'
#' @param path_to_mtx a file path to the .mtx file
#' @param dest_folder the folder in which to store the FBMs (default -- directory of .mtx file).
#'
#' @return the full file paths to the backing file (cell-gene matrix) and transposed backing file (gene-cell matrix).
#' @export
#'
#' @examples
#' path_to_mtx <- "/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/raw/GSE120861_at_scale_screen.exprs.mtx"
#' dest_folder <- "/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed"
create_fbm_from_mtx <- function(path_to_mtx, dest_folder = NULL) {
  if (is.null(dest_folder)) {
    backingfile_t <- gsub(".mtx", "_t", path_to_mtx)
    backingfile <- gsub(".mtx", "", path_to_mtx)
  } else {
    file_name <- tail(strsplit(path_to_mtx, "/")[[1]], 1)
    file_name <- gsub(".mtx", "", file_name)
    backingfile_t <- paste0(dest_folder, "/", file_name, "_t")
    backingfile <- paste0(dest_folder, "/", file_name)
  }

  header <- read.table(path_to_mtx, nrows = 1, skip = 1, header = FALSE, fill = TRUE)
  n_genes <- header[[1]]
  n_cells <- header[[2]]

  cat("Creating a file-backed matrix on disk. \n")
  expression_FBM_t <- FBM(nrow = n_genes, ncol = n_cells, type = "unsigned short",
                                 backingfile = backingfile_t,
                                 create_bk = TRUE, is_read_only = FALSE, init = 0)
  expression_FBM_t_metadata <- list(nrow = n_genes, ncol = n_cells, type = "unsigned short", backingfile = backingfile_t)
  # Define the call-back function and run the read-and-write.
  write_to_FBM <- function(chunk, pos) {
    cells_in_chunk <- chunk %>% pull(cell) %>% unique() %>% sort()
    cat(paste0("Processing cells ", cells_in_chunk[1], "-", cells_in_chunk[length(cells_in_chunk)], " of ", n_cells, ".\n"))
    for (cell in cells_in_chunk) {
      expressions_per_cell <- chunk %>% filter(cell == !!cell)
      expression_FBM_t[expressions_per_cell %>% pull(gene), cell] <- expressions_per_cell %>% pull(expression)
    }
  }

  read_delim_chunked(file = path_to_mtx, delim = " ",
                     callback = SideEffectChunkCallback$new(write_to_FBM),
                     skip = 2,
                     chunk_size = 1000000,
                     col_names = c("gene", "cell", "expression"),
                     col_types = "iii",
                     progress = FALSE)
  # Finally, save the transpose of the created matrix
  cat("Calculating and saving the transpose of the gene-cell expression matrix. \n")
  big_transpose(expression_FBM_t, backingfile = backingfile)
  expression_FBM_metadata <- list(nrow = n_cells, ncol = n_genes, type = "unsigned short", backingfile = backingfile)
  return(list(expression_FBM_metadata = expression_FBM_metadata, expression_FBM_t_metadata = expression_FBM_t_metadata))
}


#' Load FBM
#'
#' @param fbm_metadata a "metadata" list returned by the create_fbm_from_mtx function.
#'
#' @return an FBM
#' @export
#'
#' @examples
#' fbm_metadata <- readRDS("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/expression_FBM_metadata.rds")
#' load_fbm(fbm_metadata)
load_fbm <- function(fbm_metadata) {
  out <- FBM(nrow = fbm_metadata$nrow, ncol = fbm_metadata$ncol, type = fbm_metadata$type, backingfile = fbm_metadata$backingfile, is_read_only = TRUE, create_bk = FALSE)
  return(out)
}


#' Create dictionary
#'
#' Create a dictionary that maps named elements (gRNAs, genes, or gRNA pairs) to files for a given "pod" size.
#'
#' @param ids the names of elements (gRNAs or genes)
#' @param pod_size size of the pods
#' @param precomp_base_filename base name of the file (specified with full file path) in which to store the precomputation.
#'
#' @return
create_dictionary <- function(ids, pod_size) {
  n_pods_minus_1 <- floor(length(ids)/pod_size)
  pod_id <- rep(1:n_pods_minus_1, each = pod_size)
  pod_id_append <- rep(n_pods_minus_1 + 1, times = length(ids) - length(pod_id))
  pod_id <- c(pod_id, pod_id_append)
  out <- tibble(id = ids, pod_id = pod_id)
}


#' Create and store dictionaries
#'
#' We require a few bookkeeping files (that we call "dictionaries") for the gRNA and gene precomputation steps. This file creates those dictionaries and stores them in the appropriate locations on disk.
#'
#' @param gRNA_gene_pairs a tible containing a the gene-gRNA pairs to analyze using sceptre. The column names should be "gene_id" and "gRNA_id."
#' @param gene_precomp_dir the directory in which to store the gene precomputations
#' @param gRNA_precomp_dir the directory in which to store the gRNA precomputations
#' @param results_dir the directory in which to store the results
#' @param pod_sizes an integer vector with three named elements: gRNA, gene, and pair. These elements give the sizes of the respective "pods."
#'
#' @return a character vector containing the file paths to the gene, gRNA, and results dictionaries.
#' @export
create_and_store_dictionaries <- function(gRNA_gene_pairs, gene_precomp_dir, gRNA_precomp_dir, results_dir, pod_sizes) {
  # 1. genes
  genes <- unique(gRNA_gene_pairs$gene_id)
  gene_dictionary <- create_dictionary(ids = genes, pod_size = pod_sizes[["gene"]]) %>% mutate(offset_file = paste0(gene_precomp_dir, "/gene_offsets_", pod_id, ".fst") %>% factor(), dispersion_file = paste0(gene_precomp_dir, "/gene_dispersion_", pod_id, ".rds") %>% factor())
  gene_dictionary_fp <- paste0(gene_precomp_dir, "/gene_dictionary.fst")
  write.fst(gene_dictionary, gene_dictionary_fp)

  # 2. gRNA
  gRNAs <- unique(gRNA_gene_pairs$gRNA_id)
  gRNA_dictionary <- create_dictionary(ids = gRNAs, pod_size = pod_sizes[["gRNA"]]) %>% mutate(precomp_file = factor(paste0(gRNA_precomp_dir, "/gRNA_precomp_", pod_id, ".fst")))
  gRNA_dictionary_fp <- paste0(gRNA_precomp_dir, "/gRNA_dictionary.fst")
  write.fst(gRNA_dictionary, gRNA_dictionary_fp)

  # 3. gRNA-gene pairs
  pairs_dictionary <- create_dictionary(ids = gRNA_gene_pairs$gRNA_id, pod_size = pod_sizes[["pair"]]) %>% rename(gRNA_id = id) %>% mutate(gene_id = gRNA_gene_pairs$gene_id) %>% select(gRNA_id, gene_id, pod_id) %>% mutate(result_file = paste0(results_dir, "/result_", pod_id, ".fst") %>% factor())
  pairs_dictionary_fp <- paste0(results_dir, "/results_dictionary.fst")
  write.fst(pairs_dictionary, pairs_dictionary_fp)

  out <- list(fps = c(gene = gene_dictionary_fp, gRNA = gRNA_dictionary_fp, pairs = pairs_dictionary_fp),
              n_pods = c(gene = gene_dictionary$pod_id[nrow(gene_dictionary)], gRNA = gRNA_dictionary$pod_id[nrow(gRNA_dictionary)], pairs = pairs_dictionary$pod_id[nrow(pairs_dictionary)]))
  return(out)
}


#' Run gRNA precomputation at scale
#'
#' This function runs the gRNA precomputation on a selected "pod" of gRNAs (as identified by the pod_id). It stores the result in the gRNA precomputation directory.
#'
#' @param pod_id ID of the pod for which to do the precomputation
#' @param gRNA_precomp_dir file path to the gRNA precomputation directory
#' @param gRNA_indicator_matrix_fp filepath to the gRNA indicator matrix
#' @param covariate_matrix the cell-specific covariate matrix
#' @param cell_subset an integer vector specifying the cells on which to fit the model
#'
#' @return NULL
#' @export
#'
#' @examples
#' pod_id <- 1
#' gRNA_precomp_dir <- "/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/precomp/gRNA"
#' gRNA_indicator_matrix_fp <- "/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/gRNA_indicators.fst"
#' covariate_matrix <- read.fst("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/cell_covariate_model_matrix.fst")
#' cell_subset <- readRDS("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/cells_to_keep.rds")
run_gRNA_precomputation_at_scale <- function(pod_id, gRNA_precomp_dir, gRNA_indicator_matrix_fp, covariate_matrix, cell_subset = NULL, log_dir = NULL) {
  if (!is.null(log_dir)) activate_sink(paste0(log_dir, "/gRNA_precomp_", pod_id, ".Rout"))
  # subset covariate matrix according to cell subset
  if (!is.null(cell_subset)) covariate_matrix <- covariate_matrix[cell_subset,]
  # determine the gRNAs on which to run the precomputation
  gRNA_dictionary <- read.fst(paste0(gRNA_precomp_dir, "/gRNA_dictionary.fst")) %>% filter(pod_id == !!pod_id)
  gRNA_ids <- gRNA_dictionary %>% pull(id)
  # run the precomputation for each of these gRNAs
  out <- sapply(gRNA_ids, function(gRNA_id) {
    cat(paste0("Running precomputation for gRNA ", gRNA_id, ".\n"))
    gRNA_indicators <- read.fst(path = gRNA_indicator_matrix_fp, columns = gRNA_id) %>% pull() %>% as.integer()
    if (!is.null(cell_subset)) gRNA_indicators <- gRNA_indicators[cell_subset]
    run_gRNA_precomputation(gRNA_indicators, covariate_matrix)
  }) %>% as_tibble()
  # save the result
  precomp_matrix_fp <- (gRNA_dictionary %>% pull(precomp_file))[1] %>% as.character
  write.fst(out, precomp_matrix_fp)
  if (!is.null(log_dir)) deactivate_sink()
}


#' Run gene precomputation at scale
#'
#' @param pod_id ID of the gene pod on which to run the precomputation
#' @param gene_precomp_dir file path to the gene precomputation directory
#' @param cell_gene_expression_matrix an FBM containing the expression data (rows cells, columns genes)
#' @param ordered_gene_ids the gene IDs (i.e., names)
#' @param covariate_matrix the cell-specific covariate matrix
#' @param cell_subset (optional) integer vector identifying the cells to use in the model
#' @param gene_dispersions (optional) a vector of already-estimated (or known) gene dispersions
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
run_gene_precomputation_at_scale <- function(pod_id, gene_precomp_dir, cell_gene_expression_matrix, ordered_gene_ids, covariate_matrix, cell_subset = NULL, gene_dispersions = NULL, log_dir = NULL) {
  if (!is.null(log_dir)) activate_sink(paste0(log_dir, "/gene_precomp_", pod_id, ".Rout"))
  # subset covariate matrix by rows
  if (!is.null(cell_subset)) covariate_matrix <- covariate_matrix[cell_subset,]
  # obtain the genes on which to preform the precomputation
  gene_dictionary <- read.fst(paste0(gene_precomp_dir, "/gene_dictionary.fst")) %>% filter(pod_id == !!pod_id)
  gene_ids <- gene_dictionary %>% pull(id)
  integer_ids <- sapply(gene_ids, function(i) which(i == ordered_gene_ids)) %>% as.integer()

  # Run the precomputations
  precomps <- map(1:length(gene_ids), function(i) {
    cat(paste0("Running precomputation for gene ", gene_ids[i] ,".\n"))
    integer_id <- integer_ids[i]
    expressions <- cell_gene_expression_matrix[,integer_id]
    if (!is.null(cell_subset)) expressions <- expressions[cell_subset]
    run_gene_precomputation(expressions = expressions, covariate_matrix = covariate_matrix, gene_precomp_dispersion = gene_dispersions[i])
  })
  names(precomps) <- gene_ids
  out_offsets <- map_dfc(precomps, function(l) l$gene_precomp_offsets)
  out_dispersions <- map_dbl(precomps, function(l) l$gene_precomp_dispersion)

  # save the precomputations
  offset_save_fp <- (gene_dictionary %>% pull(offset_file))[1] %>% as.character()
  dispersion_save_fp <- (gene_dictionary %>% pull(dispersion_file))[1] %>% as.character()
  write.fst(out_offsets, offset_save_fp)
  saveRDS(out_dispersions, dispersion_save_fp)
  if (!is.null(log_dir)) deactivate_sink()
}


#' Run gRNA-gene pair analysis at scale
#'
#' Runs the gene-gRRNA pair ananysis across an entire pod of gene-gRNA pairs.
#'
#' @param pod_id id of the pod to compute
#' @param gene_precomp_dir the directory containing the results of the gene precomputation
#' @param gRNA_precomp_dir the directory containing the results of the gRNA precomputation
#' @param results_dir directory in which to store the results
#' @param cell_gene_expression_matrix an FBM containing the cell-by-gene expressions
#' @param ordered_gene_ids the gene ids of the cell-gene expression matrix
#' @param gRNA_indicator_matrix_fp a file-path to the gRNA indicator matrix (assumed to be stored as a .fst file)
#' @param covariate_matrix the cell-covariate matrix
#' @param cell_subset (optional) an integer vector containing the cells to subset
#' @param B number of bootstrap resamples (default 500)
#' @param log_dir (optional) directory in which to sink the log file
#' @param seed (optional) seed to pass to the randomization algorithm
#'
#' @return NULL
#' @export
#'
#' @examples
run_gRNA_gene_pair_analysis_at_scale <- function(pod_id, gene_precomp_dir, gRNA_precomp_dir, results_dir, cell_gene_expression_matrix, ordered_gene_ids, gRNA_indicator_matrix_fp, covariate_matrix, cell_subset = NULL, B = 500, log_dir = NULL, seed = NULL) {
  if (!is.null(log_dir)) activate_sink(paste0(log_dir, "/result_", pod_id, ".Rout"))

  results_dict <- read.fst(paste0(results_dir, "/results_dictionary.fst")) %>% filter(pod_id == !!pod_id)
  gene_dict <- read.fst(paste0(gene_precomp_dir, "/gene_dictionary.fst"))
  gRNA_dict <- read.fst(paste0(gRNA_precomp_dir, "/gRNA_dictionary.fst"))

  p_vals <- sapply(1:nrow(results_dict), function(i) {
    curr_gene <- results_dict[[i, "gene_id"]]
    curr_gRNA <- results_dict[[i, "gRNA_id"]]
    cat(paste0("Running distilled CRT on gene ", curr_gene, " and gRNA ", curr_gRNA, ".\n"))

    # Determine the file locations
    gene_precomp_locs <- filter(gene_dict, id == curr_gene)
    gene_offset_loc <- gene_precomp_locs %>% pull(offset_file) %>% as.character
    gene_dispersion_loc <- gene_precomp_locs %>% pull(dispersion_file) %>% as.character
    gRNA_prcomp_loc <- filter(gRNA_dict, id == curr_gRNA) %>% pull(precomp_file) %>% as.character()

    # Load the appropriate data from disk into memory
    gene_precomp_offsets <- read.fst(path = gene_offset_loc, columns = curr_gene) %>% pull()
    gene_precomp_dispersion <- readRDS(file = gene_dispersion_loc)[[curr_gene]]
    gRNA_precomp <- read.fst(path = gRNA_prcomp_loc, columns = curr_gRNA) %>% pull()
    expressions <- cell_gene_expression_matrix[, which(curr_gene == ordered_gene_ids)]
    gRNA_indicators <- read.fst(path = gRNA_indicator_matrix_fp, columns = curr_gRNA) %>% pull() %>% as.integer()

    # subset by cell id if necessary
    if (!is.null(cell_subset)) {
      expressions <- expressions[cell_subset]
      gRNA_indicators <- gRNA_indicators[cell_subset]
    }

    # Run the dCRT
    run_sceptre_using_precomp(expressions, gRNA_indicators, gRNA_precomp, gene_precomp_dispersion, gene_precomp_offsets, B, seed)
  })

  # Create and save the result dataframe
  out <- results_dict %>% summarize(gene_id = gene_id, gRNA_id = gRNA_id) %>% mutate(p_value = p_vals)
  out_fp <- (results_dict %>% pull(result_file))[1] %>% as.character()
  write.fst(out, out_fp)
  if (!is.null(log_dir)) deactivate_sink()
}


collect_results <- function(result_dir) {
  # map-reduce to obtain a single result file
  # perform BH to return the adjusted p-values as well.
}

# A convenience function to run the entire sceptre algorithm in R
run_sceptre_at_scale_simple <- function() {

}

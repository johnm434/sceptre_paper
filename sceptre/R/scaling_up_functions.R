# scaling up package
# This script provides functions to scale up sceptre to large datasets.

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
#' @export
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
  gene_dictionary <- create_dictionary(ids = genes, pod_size = pod_sizes[["gene"]]) %>% mutate(offset_file = paste0(gene_precomp_dir, "/gene_offsets_", pod_id, ".fst") %>% factor(), dispersion_file = paste0(gene_precomp_dir, "/gene_dispersion_", pod_id, ".fst") %>% factor())
  gene_dictionary_fp <- paste0(gene_precomp_dir, "/gene_dictionary.fst")
  write.fst(gene_dictionary, gene_dictionary_fp)

  # 2. gRNA
  gRNAs <- unique(gRNA_gene_pairs$gRNA_id)
  gRNA_dictionary <- create_dictionary(ids = gRNAs, pod_size = pod_sizes[["gRNA"]]) %>% mutate(precomp_file = factor(paste0(gRNA_precomp_dir, "/gRNA_precomp_", pod_id, ".fst")))
  gRNA_dictionary_fp <- paste0(gRNA_precomp_dir, "/gRNA_dictionary.fst")
  write.fst(gRNA_dictionary, gRNA_dictionary_fp)

  # 3. gRNA-gene pairs
  pairs_dictionary <- create_dictionary(ids = gRNA_gene_pairs$gRNA_id, pod_size = pod_sizes[["pair"]]) %>% rename(gRNA_id = id) %>% mutate(gene_id = gRNA_gene_pairs$gene_id) %>% select(gRNA_id, gene_id, pod_id) %>% mutate(result_file = paste0(results_dir, "/result_", pod_id, ".fst"))
  pairs_dictionary_fp <- paste0(results_dir, "/results_dictionary.fst")
  write.fst(pairs_dictionary, pairs_dictionary_fp)

  out <- list(fps = c(gene = gene_dictionary_fp, gRNA = gRNA_dictionary_fp, pairs = pairs_dictionary_fp),
              n_pods = c(gene = gene_dictionary$pod_id[nrow(gene_dictionary)], gRNA = gRNA_dictionary$pod_id[nrow(gRNA_dictionary)], pairs = pairs_dictionary$pod_id[nrow(pairs_dictionary)]))
  return(out)
}


#' Run gRNA precomputation at scale
#'
#' @param gRNA_ids the ids of gRNAs on which to compute the precomputation
#' @param gRNA_indicator_matrix_fp file path the the gRNA indicator matrix (assumed to be stored in .fst format)
#' @param covariate_matrix the covariate matrix
#' @param precomp_matrix_fp the name of the file (including the file path) in which to save the precomputation result
#' @param cell_subset (optional) an integer vector containing the cells to use in the precomputation
#'
#' @return NULL
#' @export
#'
#' @examples
#' gRNA_dict <- read.fst("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/precomp/gRNA/gRNA_dictionary.fst")
#' gRNA_ids <- filter(gRNA_dict, pod_id == 1) %>% pull(id) %>% head(20)
#' precomp_matrix_fp <- as.character((filter(gRNA_dict, pod_id == 1) %>% pull(file_name))[1])
#' gRNA_indicator_matrix_fp <- "/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/gRNA_indicators.fst"
#' covariate_matrix <- read.fst("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/cell_covariate_model_matrix.fst")
#' run_gRNA_precomputation_at_scale(gRNA_ids[1:20], gRNA_indicator_matrix_fp, covariate_matrix, precomp_matrix_fp)

run_gRNA_precomputation_at_scale <- function(gRNA_ids, gRNA_indicator_matrix_fp, covariate_matrix, precomp_matrix_fp, cell_subset = NULL) {
  if (!is.null(cell_subset)) covariate_matrix <- covariate_matrix[cell_subset,]
  # run the precomputation on the gRNAs
  out <- sapply(gRNA_ids, function(gRNA_id) {
    cat(paste0("Running precomputation for gRNA ", gRNA_id, ".\n"))
    gRNA_indicators <- read.fst(path = gRNA_indicator_matrix_fp, columns = gRNA_id) %>% pull() %>% as.integer()
    if (!is.null(cell_subset)) gRNA_indicators <- gRNA_indicators[cell_subset]
    run_gRNA_precomputation(gRNA_indicators, covariate_matrix)
  }) %>% as_tibble()
  # save the result
  write.fst(out, precomp_matrix_fp)
  return(NULL)
}

# gene_dict <- read.fst("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/precomp/gene/gene_dictionary.fst")
# genes_to_precompute <- filter(gene_dict, pod_id == 1) %>% pull(id) %>% head(5)
# offset_save_fp <- "/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/precomp/gene/offset_precomp_1.fst"
# dispersion_save_fp <- "/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/precomp/gene/dispersion_precomp_1.rds"
# covariate_matrix <- read.fst("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/cell_covariate_model_matrix.fst")
# exp_matrix_info <- readRDS("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/expression_FBM_metadata.rds")
# exp_matrix <- FBM(nrow = exp_matrix_info$nrow, ncol = exp_matrix_info$ncol, type = exp_matrix_info$type, create_bk = FALSE, is_read_only = TRUE, backingfile = exp_matrix_info$backingfile)
# all_ordered_gene_ids <- readRDS("/Volumes/tims_new_drive/research/sceptre_files/data/gasperini/processed/ordered_gene_ids.RDS")

run_gene_precomputation_at_scale <- function(genes_to_precompute, exp_matrix, all_ordered_gene_ids, covariate_matrix, offset_save_fp, dispersion_save_fp, gene_dispersions = NULL, cell_subset = NULL) {
  if (!is.null(cell_subset)) covariate_matrix <- covariate_matrix[cell_subset,]
  integer_ids <- sapply(genes_to_precompute, function(i) which(i == all_ordered_gene_ids)) %>% as.integer()
  precomps <- map(1:length(genes_to_precompute), function(i) {
    cat(paste0("Running precomputation for gene ", genes_to_precompute[i] ,".\n"))
    integer_id <- integer_ids[i]
    expressions <- exp_matrix[,integer_id]
    if (!is.null(cell_subset)) expressions <- expressions[cell_subset]
    run_gene_precomputation(expressions = expressions, covariate_matrix = covariate_matrix, gene_precomp_dispersion = gene_dispersions[i], NULL)
  })
  names(precomps) <- genes_to_precompute
  out_offsets <- map_dfc(precomps, function(l) l$gene_precomp_offsets)
  out_dispersions <- map_dbl(precomps, function(l) l$gene_precomp_dispersion)
  write.fst(out_offsets, offset_save_fp)
  saveRDS(out_dispersions, dispersion_save_fp)
}

# Pre-process data
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


#' Extract column from CSC (compreseed sparse column) representation of matrix
#'
#' @param column_no integer vector; the columns to extract from the matrix
#' @param dat the vector of UMI counts
#' @param indices the "external" indices corresponding to row ids
#' @param ind_ptr the "internal" indices used to locate the external indices and dat for a given column
#' @param n_rows number of rows of matrix
#' @param zero_based_idx (boolean) do the indices and ind_ptr vectors use zero-based indexing?
#'
#' @return the (dense) columns of the matrix corresponding to columns column_no
#' @export
extract_column_from_csc <- function(column_no, dat, indices, ind_ptr, n_rows, zero_based_idx = TRUE) {
  if (zero_based_idx) { # make the vectors use 1-based indexing
    indices <- indices + 1
    ind_ptr <- ind_ptr + 1
  }
  sapply(X = column_no, FUN = function(i) {
    internal_idx_range <- ind_ptr[i:(i + 1)] + c(0,-1)
    non_zero_out <- dat[internal_idx_range[1]:internal_idx_range[2]]
    external_idx <- indices[internal_idx_range[1]:internal_idx_range[2]]
    full_out <- rep(0, n_rows)
    full_out[external_idx] <- non_zero_out
    return(full_out)
  })
}

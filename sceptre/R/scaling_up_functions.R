# scaling up package
# This script provides functions to scale up sceptre to large datasets.


#' Create fbm from mtx
#'
#' @param path_to_mtx
#' @param dest_folder
#'
#' @return
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
  expression_FBM_t <- FBM(n_genes, n_cells, type = "unsigned short",
                                 backingfile = backingfile_t,
                                 create_bk = TRUE, is_read_only = FALSE, init = 0)
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
  big_transpose(expression_FBM_t, backingfile = backingfile)
}

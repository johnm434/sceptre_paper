#' Get scaled data
#'
#' Obtain the scaled data matrix from a normalized/scaled Seurat object.
#'
#' This function was copied and pased from the scMAGeCK package; Yang et al. wrote the function, not Katsevich et al.#' Single gene matrix regression
#'
#' Obtain the X and Y matrices for a regression.
#'
#' This function was copied and pased from the scMAGeCK package; Yang et al. wrote the function, not Katsevich et al.
#'
#' @param targetobj a scaled Seurat object
#' @param ngctrlgene the name of a gene or genes that serves as a negative control
#' @param indmatrix the matrix of gRNA-gene indicators
#' @param high_gene_frac remove expressions above this level
#' @param selected_genes_list generate Y matrix for these genes
#'
#' @return and X and Y matrix for regression
#' @export
single_gene_matrix_regression <- function(targetobj, ngctrlgene = c("NonTargetingControlGuideForHuman"),
                                          indmatrix = NULL, high_gene_frac = 0.01, selected_genes_list = NULL) {
  # return X matrix and Y matrix for regression note that all the ngctrlgene are merged into one
  # column, 'NegCtrl' if indmatrix is provided, the Xmat will be constructed from indmatrix
  outlier_threshold = 0.95
  rawf = targetobj@assays$RNA@counts %>% as.matrix()
  select_genes = rownames(rawf)[which(rowSums(as.matrix(rawf) != 0) >= ncol(rawf) * high_gene_frac)]
  if (is.null(selected_genes_list) == FALSE) {
    select_genes = select_genes[select_genes %in% selected_genes_list]
    if (length(select_genes) == 0) {
      stop("No genes left for regression. Check your selected gene list.")
    }
  }
  message(paste("Selected genes:", length(select_genes)))

  scalef = targetobj[["RNA"]]@data %>% as.matrix()

  if (is.null(indmatrix)) {
    select_cells = rownames(targetobj@meta.data)[which(!is.na(targetobj@meta.data$geneID))]
  } else {
    select_cells = rownames(indmatrix)
    select_cells = select_cells[select_cells %in% colnames(scalef)]
  }
  YmatT = scalef[select_genes, select_cells]

  Ymat = as.matrix(t(YmatT))  # (cells * expressed genes)
  if (is.null(indmatrix)) {
    tgf = targetobj@meta.data[select_cells, "geneID"]
    tgf[tgf %in% ngctrlgene] = "NegCtrl"
    tgphenotype = as.factor(tgf)
    Xmat = matrix(rep(0, length(select_cells) * length(unique(tgphenotype))), nrow = length(select_cells))
    rownames(Xmat) = select_cells
    colnames(Xmat) = levels(tgphenotype)
    Xmat[as.matrix(cbind(1:nrow(Xmat), as.numeric(tgphenotype)))] = 1
    Xmat[, "NegCtrl"] = 1  # set up base line
  } else {
    tgf = colnames(indmatrix)
    tgf[tgf %in% ngctrlgene] = "NegCtrl"
    tgphenotype = as.factor(tgf)

    Xmat = matrix(rep(0, length(select_cells) * length(unique(tgphenotype))), nrow = length(select_cells))
    rownames(Xmat) = select_cells
    colnames(Xmat) = levels(tgphenotype)
    for (cnl in colnames(indmatrix)) {
      cellns = which(indmatrix[, cnl] == TRUE)  #make sure indmatrix
      if (cnl %in% ngctrlgene) {
        Xmat[cellns, "NegCtrl"] = 1
      } else {
        Xmat[cellns, cnl] = 1
      }
    }
    Xmat[, "NegCtrl"] = 1
  }

  # remove outliers
  Ymat_outlier = apply(Ymat, 2, function(X) {
    return(quantile(X, probs = outlier_threshold))
  })
  outlier_mat = t(matrix(rep(Ymat_outlier, nrow(Ymat)), ncol = nrow(Ymat)))
  Ymat_corrected = ifelse(Ymat > outlier_mat, outlier_mat, Ymat)
  Ymat = Ymat_corrected

  return(list(Xmat, Ymat))
}


#' getsolvedmatrix_with_permutation_cell_label
#'
#' This function was copied and pased from the scMAGeCK package; Yang et al. wrote the function, not Katsevich et al.
#'
#' @param Xm the X matrix for the regression
#' @param Ym the Y matrix
#' @param lambda the ridge penalty
#' @param npermutation number of permutations to run
#'
#' @return
#' @export
getsolvedmatrix_with_permutation_cell_label <- function(Xm, Ym, lambda = 0.01, npermutation = 1000) {
  Amat_ret = getsolvedmatrix(Xm, Ym, lambda = lambda)
  Amat_ret_higher = matrix(rep(0, ncol(Amat_ret) * nrow(Amat_ret)), nrow = nrow(Amat_ret))
  rownames(Amat_ret_higher) = rownames(Amat_ret)
  colnames(Amat_ret_higher) = colnames(Amat_ret)
  # permute N times randomly shuffle cell labels
  for (npm in 1:npermutation) {
    if (npm%%100 == 0) {
      message(paste("Permutation:", npm, "/", npermutation, "..."))
    }
    cells_shu = sample(rownames(Ym), nrow(Ym))
    Xm_s = Xm[cells_shu, ]
    Ym_s = Ym  # [cells_shu,]
    rownames(Ym_s) = cells_shu
    Amat_random = getsolvedmatrix(Xm_s, Ym_s, lambda = lambda)

    Amat_ret_higher = Amat_ret_higher + (abs(Amat_random) > abs(Amat_ret)) * 1
    # browser()
  }
  Amat_ret_higher = Amat_ret_higher/npermutation
  return(list(Amat_ret, Amat_ret_higher))
}


#' Get solved matrix
#'
#' This function was copied and pased from the scMAGeCK package; Yang et al. wrote the function, not Katsevich et al.
#'
#' @param Xm the X matrix for the regression
#' @param Ym the Y matrix
#' @param lambda the ridge penalty
#'
#' @return
#' @export
getsolvedmatrix <- function(Xm, Ym, lambda = 0.01) {
  # Amat=solve(Xmat,Ymat) # solve AX=B, or Xmat * A =Ymat
  TMmat_g = (t(Xm) %*% Xm) + lambda * diag(ncol(Xm))

  Amat_g = solve(TMmat_g) %*% t(Xm) %*% Ym
  return(Amat_g)
}


#' Run scMaGECK analysis
#'
#' @param expression_matrix a matrix of cell-by-gene expressions
#' @param gRNA_indic_matrix a matrix of cell-by-gRNA indicators
#'
#' @return p-values for each gRNA-gene pair
#' @export
#'
#' @examples
#' r <- simulate_crispr_screen_data(n_cells = 1000,
#' grna_mean_prob = 0.2,
#' mRNA_mean_expression = 40,
#' covariate_effects_gRNA = rep(4, 10),
#' covariate_effects_gene = rep(2, 50),
#' zero_inflation = 0,
#' neg_binom_size = 2
#' )
#' expression_matrix <- r$cell_by_gene_expressions
#' gRNA_indic_matrix <- r$cell_by_enhancer_perturbation_indicators
#' run_scMaGECK_analysis(expression_matrix, gRNA_indic_matrix)
run_scMaGECK_analysis <- function(expression_matrix, gRNA_indic_matrix) {
  expression_matrix_t <- t(expression_matrix)
  row.names(expression_matrix_t) <- paste0("gene-", 1:nrow(expression_matrix_t))
  colnames(expression_matrix_t) <- paste0("cell-", 1:ncol(expression_matrix_t))
  row.names(gRNA_indic_matrix) <- paste0("cell-", 1:nrow(gRNA_indic_matrix))
  colnames(gRNA_indic_matrix) <- paste0("gRNA-", 1:ncol(gRNA_indic_matrix))
  s <- CreateSeuratObject(expression_matrix_t)

  targetobj <- NormalizeData(s)
  GENE_FRAC <- 0.00
  SELECT_GENE <- row.names(expression_matrix_t)
  ngctrlgenelist <- colnames(gRNA_indic_matrix)[ncol(gRNA_indic_matrix)]
  ind_matrix <- gRNA_indic_matrix

  mat_for_single_reg = single_gene_matrix_regression(targetobj,
                                                     selected_genes_list = SELECT_GENE, ngctrlgene = ngctrlgenelist,
                                                     indmatrix = ind_matrix, high_gene_frac = GENE_FRAC)
  Xmat = mat_for_single_reg[[1]]
  Ymat = mat_for_single_reg[[2]]
  LAMBDA <- 0.01
  Amat_pm_lst = getsolvedmatrix_with_permutation_cell_label(Xmat, Ymat, lambda = LAMBDA, npermutation = 500)
  Amat = Amat_pm_lst[[1]]
  Amat_pval = Amat_pm_lst[[2]]
  return(Amat_pval)
}

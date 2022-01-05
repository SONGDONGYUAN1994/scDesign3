#' Construct the input data (covaraite matrix and expression matrix)
#'
#' This function constructs the input data for \code{\link{fit_marginal}}.
#'
#' This function takes a \code{SingleCellExperiment} object as the input.
#' Based on users' choice, it constructs the matrix of covaraites
#' (explainary variables) and the expression matrix (e.g., count matrix for scRNA-seq).
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param assay_use A string which indicates the assay you will use in the sce. Default is 'counts'.
#' @param covariate_use A string of the primary covariate.
#' Must be one of 'celltype', 'pseudotime' or 'spatial'.
#' @param celltype A string of the name of cell type variable in the \code{colData} of the sce. Default is 'cee_type'.
#' @param pseudotime A string or a string vector of the name of pseudotime and (if exist)
#' multiple lineages. Default is NULL.
#' @param spatial A length two string vector of the names of spatial coordinates. Defualt is NULL.
#' @param other_covariates A string or a string vector of the other covaraites you want to include in the data.
#' @param ncell The number of cell you want to simulate. Default is \code{dim(sce)[2]} (the same number as the input data).
#' If an arbitrary number is provided, the fucntion will use Vine Copula to simulate a new covaraite matrix.
#'
#' @return A list with the components:
#' \describe{
#'   \item{\code{count_mat}}{The expression matrix}
#'   \item{\code{dat}}{The original covariate matrix}
#'   \item{\code{newCovariate}}{The simulated new covariate matrix}
#' }
#'
#' @export construct_data
construct_data <- function(sce,
                          assay_use = "counts",
                          covariate_use,
                          celltype,
                          pseudotime,
                          spatial,
                          other_covariates,
                          ncell = dim(sce)[2]) {
  ## Extract expression matrix
  count_mat <-
    t(as.matrix(SummarizedExperiment::assay(sce, assay_use)))

  ## Extract col data
  coldata_mat <- data.frame(SummarizedExperiment::colData(sce))



  ## Extract pseudotime / cell type / spatial
  if (!is.null(celltype)) {
    celltype <- as.matrix(coldata_mat[, celltype, drop = FALSE])
  }

  if (!is.null(pseudotime)) {
    pseudotime <- as.matrix(coldata_mat[, pseudotime, drop = FALSE])
  }

  if (!is.null(spatial)) {
    spatial <- as.matrix(coldata_mat[, spatial, drop = FALSE])
  }


  if (covariate_use == "celltype") {
    dat <- data.frame(celltype)
    dat$cell_type <- as.factor(dat$cell_type)
  } else if (covariate_use == "pseudotime") {
    n_l <- dim(pseudotime)[2]
    dat <- data.frame(pseudotime)
  } else if (covariate_use == "spatial") {
    dat <- data.frame(spatial)
  } else {
    stop("Covairate_use must be one of 'celltype', 'pseudotime' or 'spatial'!")
  }

  ## Convert NA to -1
  #pseudotime[is.na(pseudotime)] <- -1

  ## dat is the input covariate matrix
  if (!is.null(other_covariates)) {
    other_covariates <- as.matrix(coldata_mat[, other_covariates, drop = FALSE])
    dat <- cbind(dat, other_covariates)
    if("condition" %in% colnames(other_covariates)){
      dat$condition <- as.factor(dat$condition)
    }
  }

  ## check if user wants to simulate new number of cells
  if(ncell != dim(dat)[1]){
    newCovariate <- as.data.frame(simuCovariateMat(dat,ncell))
  }else{
    newCovariate <- NULL
  }

  return(list(count_mat = count_mat, dat = dat, newCovariate = newCovariate))
}


## Simulate covariate matrix
simuCovariateMat <- function(covariate_mat,
                             n_cell_new = 50000) {

  n_cell_ori <- dim(covariate_mat)[1]
  n_covraite_ori <- dim(covariate_mat)[2]

  if_factor_exist <- sum(sapply(covariate_mat, is.factor))
  if_numeric_exist <- sum(sapply(covariate_mat, is.numeric))



  df <- covariate_mat

  if(if_factor_exist) {
    df_all <- dplyr::mutate(df, discrete_group = interaction(dplyr::mutate_if(is.factor), sep = "-"))

    df_list <- split(df_all, df_all$discrete_group)

    group_prop <-  table(df_all$discrete_group)/dim(df_all)[1]
    group_name <- names(group_prop)
    group_n_new <- stats::rmultinom(1, size = n_cell_new, prob = group_prop)
    if(if_numeric_exist) {
      new_dat_list <- pbmcapply::pbmcmapply(FUN = function(df, n) {
        df <- dplyr::select(df, -"discrete_group")
        df_numeric <- dplyr::select_if(df, is.numeric)
        df_factor <- dplyr::select_if(df, is.factor)
        fit_kde <- rvinecopulib::vine(df_numeric, cores = 1)
        new_dat <- rvinecopulib::rvine(n, fit_kde)

        new_dat <- as.data.frame(new_dat)
        new_dat[colnames(df_factor)] <- df_factor[1, ]
        new_dat

      }, df = df_list, n = group_n_new, SIMPLIFY = FALSE)
      covariate_new <- do.call("rbind", new_dat_list)
    }
    else {
      new_dat_list <- pbmcapply::pbmcmapply(FUN = function(df, n) {
        df <- dplyr::select(df, -"discrete_group")
        df_factor <- dplyr::select_if(df, is.factor)
        df_onerow <- as.data.frame(df_factor[1, ])
        new_dat <- as.data.frame(df_onerow[rep(1, n), ])
      }, df = df_list, n = group_n_new, SIMPLIFY = FALSE)
      covariate_new <- do.call("rbind", new_dat_list)
      colnames(covariate_new) <- colnames(covariate_mat)
    }
  }
  else {
    fit_kde <- rvinecopulib::vine(df, cores = 1)
    covariate_new <- rvinecopulib::rvine(n_cell_new, fit_kde)
  }

  rownames(covariate_new) <- paste0("Cell", seq_len(n_cell_new))
  return(covariate_new)
}


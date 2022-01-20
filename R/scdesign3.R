#' The wrapper for the whole scDesign3 pipeline
#'
#' \code{scdesign3} takes the input data, fits the model and
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param assay_use A string which indicates the assay you will use in the sce. Default is 'counts'.
#' Must be one of 'celltype', 'pseudotime' or 'spatial'.
#' @param celltype A string of the name of cell type variable in the \code{colData} of the sce. Default is 'cell_type'.
#' @param pseudotime A string or a string vector of the name of pseudotime and (if exist)
#' multiple lineages. Default is NULL.
#' @param spatial A length two string vector of the names of spatial coordinates. Defualt is NULL.
#' @param other_covariates A string or a string vector of the other covaraites you want to include in the data.
#' @param ncell The number of cell you want to simulate. Default is \code{dim(sce)[2]} (the same number as the input data).
#' @param mu_formula A string of the mu parameter formula
#' @param sigma_formula A string of the sigma parameter formula
#' @param family_use A string of the marginal distribution.
#' Must be one of 'poisson', 'nb', 'zip', 'zinb' or 'gaussian'.
#' @param n_cores An integer. The number of cores to use.
#' @param usebam A logic variable. If use \code{\link[mgcv]{bam}} for acceleration.
#' @param corr_formula A string of the correlation structure.
#' @param copula A string of the copula choice. Must be one of 'gaussian' or 'vine'. Default is 'vine'.
#' @param DT A logic variable. If TRUE, perform the distributional transformation
#' to make the discrete data 'continuous'. This is useful for discrete distributions (e.g., Poisson, NB).
#' Default is TRUE.
#' @param pseudo_obs A logic variable. If TRUE, use the empirical quantiles instead of theoretical quantiles for fitting copula.
#' Default is FALSE.
#' @param family_set A string or a string vector of the bivariate copula families. Default is c("gaussian", "indep").
#' @param return_model A logic variable. If TRUE, the marginal models and copula models will be returned. Default is FALSE.
#'
#' @return A list with the components:
#' \describe{
#'   \item{\code{new_count}}{A matrix of the new simulated count (expression) matrix.}
#'   \item{\code{new_covariate}}{A data.frame of the new covariate matrix.}
#' }
#' @export scdesign3
scdesign3 <- function(sce,
                      assay_use = "counts",
                      celltype,
                      pseudotime,
                      spatial,
                      other_covariates,
                      ncell = dim(sce)[2],
                      mu_formula,
                      sigma_formula = "1",
                      family_use = "nb",
                      n_cores = 2,
                      usebam = FALSE,
                      corr_formula,
                      copula = "vine",
                      DT = TRUE,
                      pseudo_obs = FALSE,
                      family_set = c("gauss"),
                      return_model = FALSE) {
  message("Input Data Construction Start")

  input_data <- construct_data(
    sce = sce,
    assay_use = assay_use,
    celltype = celltype,
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates,
    ncell = ncell,
    corr_by = corr_formula
  )
  message("Input Data Construction End")

  message("Start Marginal Fitting")
  marginal_res <- fit_marginal(
    mu_formula = mu_formula,
    sigma_formula = sigma_formula,
    n_cores = n_cores,
    data = input_data,
    family_use = family_use,
    usebam = usebam
  )
  message("Marginal Fitting End")

  message("Start Copula Fitting")
  copula_res <- fit_copula(
    sce = sce,
    assay_use = assay_use,
    input_data = input_data$dat,
    new_covariate = input_data$new_covariate,
    marginal_list = marginal_res,
    family_use = family_use,
    copula = copula,
    family_set = family_set,
    n_cores = n_cores
  )
  message("Copula Fitting End")

  message("Start Parameter Extraction")
  para_list <- extract_para(
    sce = sce,
    marginal_list = marginal_res,
    n_cores = n_cores,
    family_use = family_use,
    new_covariate = input_data$new_covariate
  )
  message("Parameter
Extraction End")

  message("Start Generate New Data")
  new_count <- simu_new(
    sce = sce,
    mean_mat = para_list$mean_mat,
    sigma_mat = para_list$sigma_mat,
    zero_mat = para_list$zero_mat,
    quantile_mat = copula_res$new_mvu,
    n_cores = n_cores,
    family_use = family_use,
    new_covariate = input_data$new_covariate
  )
  message("New Data Generating End")

  scdesign3_res <- list(
    new_count = new_count,
    new_covariate = input_data$new_covariate,
    model_aic = copula_res$model_aic,
    model_bic = copula_res$model_bic,
    marginal_list = if (return_model)
      marginal_res
    else
      NULL,
    corr_list = if (return_model)
      copula_res$corr_list
    else
      NULL
  )
  return(scdesign3_res)
}

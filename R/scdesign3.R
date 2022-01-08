#' The wrapper for the whole scDesign3 pipeline
#'
#' \code{scdesign3} takes the input data, fits the model and
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param assay_use A string which indicates the assay you will use in the sce. Default is 'counts'.
#' @param covariate_use A string of the primary covariate.
#' Must be one of 'celltype', 'pseudotime' or 'spatial'.
#' @param celltype A string of the name of cell type variable in the \code{colData} of the sce. Default is 'cell_type'.
#' @param pseudotime A string or a string vector of the name of pseudotime and (if exist)
#' multiple lineages. Default is NULL.
#' @param spatial A length two string vector of the names of spatial coordinates. Defualt is NULL.
#' @param other_covariates A string or a string vector of the other covaraites you want to include in the data.
#' @param ncell The number of cell you want to simulate. Default is \code{dim(sce)[2]} (the same number as the input data).
#' @param predictor Default is gene. ## Fix later
#' @param mu_formula A string of the mu parameter formula
#' @param sigma_formula A string of the sigma parameter formula
#' @param family A string of the marginal distribution.
#' Must be one of 'poisson', 'nb', 'zip', 'zinb' or 'gaussian'.
#' @param n_cores An integer. The number of cores to use.
#' @param usebam A logic variable. If use \code{\link[mgcv]{bam}} for acceleration.
#' @param cor_formula A string of the correlation structure.
#' @param copula A string of the copula choice. Must be one of 'gaussian' or 'vine'. Default is 'vine'.
#' @param DT A logic variable. If TRUE, perform the distributional transformation
#' to make the discrete data 'continuous'. This is useful for discrete distributions (e.g., Poisson, NB).
#' Default is TRUE.
#' @param pseudo_obs A logic variable. If TRUE, use the empirical quantiles instead of theoretical quantiles for fitting copula.
#' Default is FALSE.
#' @param family_set A string or a string vector of the bivariate copula families. Default is c("gaussian", "indep").
#'
#' @return A list with the components:
#' \describe{
#'   \item{\code{new_count}}{A matrix of the new simulated count (expression) matrix.}
#'   \item{\code{new_covariate}}{A data.frame of the new covariate matrix.}
#' }
#' @export scdesign3
scdesign3 <- function(sce,
                      assay_use = "counts",
                      covariate_use = "celltype",
                      celltype = "cell_type",
                      pseudotime = "pseudotime",
                      spatial = c("spatial1", "spatial2"),
                      other_covariates = c("batch", "condition"),
                      ncell = dim(sce)[2],
                      predictor = "gene",
                      mu_formula = "s(pseudotime, k = 10, bs = 'cr')",
                      sigma_formula = "1",
                      family = "nb",
                      n_cores = 2,
                      usebam = FALSE,
                      cor_formula = "cell_type",
                      copula = "vine",
                      DT = TRUE,
                      pseudo_obs = FALSE,
                      family_set = c("gaussian", "indep")
                      ){

  message("Input Data Construction Start")

  input_data <- construct_data(sce = sce,
                        assay_use = assay_use,
                        celltype = celltype,
                        pseudotime = pseudotime,
                        spatial = spatial,
                        other_covariates = other_covariates,
                        ncell = ncell,
                        corr_by = cor_formula)
  message("Input Data Construction End")

  message("Start Marginal Fitting")
  marginal_res <- fit_marginal(
    mu_formula = mu_formula,
    sigma_formula = sigma_formula,
    n_cores = n_cores,
    data = input_data,
    predictor = predictor,
    family = family,
    usebam = usebam
  )
  message("Marginal Fitting End")

  message("Start Copula Fitting")
  copula_res <- fit_copula(sce = sce,
                          assay_use = assay_use,
                          input_data = input_data$dat,
                          new_covariate = input_data$new_covariate,
                          marginal_list = marginal_res,
                          family = family,
                          copula = copula,
                          cor_formula = cor_formula,
                          family_set = family_set,
                          n_cores = n_cores)
  message("Copula Fitting End")

  message("Start Generate New Data")
  res_list <- simu_new(sce = sce,
                           marginal_list = marginal_res,
                           copula_list = copula_res,
                           n_cores = n_cores ,
                           family = family,
                           new_covariate = input_data$new_covariate)
  message("New Data Generating End")
  return(list(new_count = res_list$new_count, new_covariate = input_data$new_covariate))
}

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
#' @param other_covariates A string or a string vector of the other covariates you want to include in the data.
#' @param ncell The number of cell you want to simulate. Default is \code{dim(sce)[2]} (the same number as the input data).
#' @param mu_formula A string of the mu parameter formula
#' @param sigma_formula A string of the sigma parameter formula
#' @param family_use A string of the marginal distribution.
#' Must be one of 'poisson', 'nb', 'zip', 'zinb' or 'gaussian'.
#' @param n_cores An integer. The number of cores to use.
#' @param usebam A logic variable. If use \code{\link[mgcv]{bam}} for acceleration.
#' @param corr_formula A string of the correlation structure.
#' @param copula A string of the copula choice. Must be one of 'gaussian' or 'vine'. Default is 'gaussian'. Note that vine copula may have better modeling of high-dimensions, but can be very slow when features are >1000.
#' @param fastmvn An logical variable. If TRUE, the sampling of multivariate Gaussian is done by \code{mvnfast}, otherwise by \code{mvtnorm}. Default is FALSE. It only matters for Gaussian copula.
#' @param DT A logic variable. If TRUE, perform the distributional transformation
#' to make the discrete data 'continuous'. This is useful for discrete distributions (e.g., Poisson, NB).
#' Default is TRUE. Note that for continuous data (e.g., Gaussian), DT does not make sense and should be set as FALSE.
#' @param pseudo_obs A logic variable. If TRUE, use the empirical quantiles instead of theoretical quantiles for fitting copula.
#' Default is FALSE.
#' @param family_set A string or a string vector of the bivariate copula families. Default is c("gauss", "indep").
#' @param important_feature A string or vector which indicates whether a gene will be used in correlation estimation or not. If this is a string, then
#' this string must be "auto", which indicates that the genes will be automatically selected based on the proportion of zero expression across cells
#' for each gene. Gene with zero proportion greater than 0.8 will be excluded form gene-gene correlation estimation. If this is a vector, then this should
#' be a logical vector with length equal to the number of genes in \code{sce}. \code{TRUE} in the logical vector means the corresponding gene will be included in
#' gene-gene correlation estimation and \code{FALSE} in the logical vector means the corresponding gene will be excluded from the gene-gene correlation estimation.
#' The default value for is a vector with length equal to the number of inputted genes and every value equals to \code{TRUE}.
#' @param nonnegative A logical variable. If TRUE, values < 0 in the synthetic data will be converted to 0. Default is TRUE (since the expression matrix is nonnegative).
#' @param nonzerovar A logical variable. If TRUE, for any gene with zero variance, a cell will be replaced with 1. This is designed for avoiding potential errors, for example, PCA.
#' @param return_model A logic variable. If TRUE, the marginal models and copula models will be returned. Default is FALSE.
#' @param parallelization A string indicating the specific parallelization function to use.
#' Must be one of 'mcmapply', 'bpmapply', or 'pbmcmapply', which corresponds to the parallelization function in the package
#' \code{parallel},\code{BiocParallel}, and \code{pbmcapply} respectively. The default value is 'mcmapply'.
#' @param BPPARAM A \code{MulticoreParam} object or NULL. When the parameter parallelization = 'mcmapply' or 'pbmcmapply',
#' this parameter must be NULL. When the parameter parallelization = 'bpmapply',  this parameter must be one of the
#' \code{MulticoreParam} object offered by the package 'BiocParallel. The default value is NULL.
#' @param trace A logic variable. If TRUE, the warning/error log and runtime for gam/gamlss
#' will be returned, FALSE otherwise. Default is FALSE.
#' @return A list with the components:
#' \describe{
#'   \item{\code{new_count}}{A matrix of the new simulated count (expression) matrix.}
#'   \item{\code{new_covariate}}{A data.frame of the new covariate matrix.}
#'   \item{\code{model_aic}}{The model AIC.}
#'   \item{\code{marginal_list}}{A list of marginal regression models if return_model = TRUE.}
#'   \item{\code{corr_list}}{A list of correlation models (conditional copulas) if return_model = TRUE.}
#' }
#' @examples
#' data(example_sce)
#' my_simu <- scdesign3(
#' sce = example_sce,
#' assay_use = "counts",
#' celltype = "cell_type",
#' pseudotime = "pseudotime",
#' spatial = NULL,
#' other_covariates = NULL,
#' mu_formula = "s(pseudotime, bs = 'cr', k = 10)",
#' sigma_formula = "s(pseudotime, bs = 'cr', k = 3)",
#' family_use = c(rep("nb", 5), rep("zip", 5)),
#' n_cores = 2,
#' usebam = FALSE,
#' corr_formula = "pseudotime",
#' copula = "vine",
#' DT = TRUE,
#' pseudo_obs = FALSE,
#' ncell = 1000,
#' return_model = FALSE
#' )
#'
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
                      copula = "gaussian",
                      fastmvn = FALSE,
                      DT = TRUE,
                      pseudo_obs = FALSE,
                      family_set = c("gauss", "indep"),
                      important_feature = rep(TRUE, dim(sce)[1]),
                      nonnegative = TRUE,
                      nonzerovar = TRUE,
                      return_model = FALSE,
                      parallelization = "mcmapply",
                      BPPARAM = NULL,
                      trace = FALSE) {
  message("Input Data Construction Start")

  input_data <- construct_data(
    sce = sce,
    assay_use = assay_use,
    celltype = celltype,
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates,
    ncell = ncell,
    corr_by = corr_formula,
    parallelization = "mcmapply",
    BPPARAM = BPPARAM
  )
  message("Input Data Construction End")

  message("Start Marginal Fitting")
  marginal_res <- fit_marginal(
    mu_formula = mu_formula,
    sigma_formula = sigma_formula,
    n_cores = n_cores,
    data = input_data,
    family_use = family_use,
    usebam = usebam,
    parallelization = parallelization,
    BPPARAM = BPPARAM,
    trace = trace
  )
  message("Marginal Fitting End")

  message("Start Copula Fitting")
  copula_res <- fit_copula(
    sce = sce,
    assay_use = assay_use,
    input_data = input_data$dat,
    new_covariate = input_data$newCovariate,
    marginal_list = marginal_res,
    family_use = family_use,
    copula = copula,
    family_set = family_set,
    n_cores = n_cores,
    important_feature = important_feature,
    parallelization = parallelization,
    BPPARAM = BPPARAM
  )
  message("Copula Fitting End")

  message("Start Parameter Extraction")
  para_list <- extract_para(
    sce = sce,
    assay_use = assay_use,
    marginal_list = marginal_res,
    n_cores = n_cores,
    family_use = family_use,
    new_covariate = input_data$newCovariate,
    parallelization = parallelization,
    BPPARAM = BPPARAM,
    data = input_data$dat
  )
  message("Parameter
Extraction End")

  message("Start Generate New Data")
  new_count <- simu_new(
    sce = sce,
    assay_use= assay_use,
    mean_mat = para_list$mean_mat,
    sigma_mat = para_list$sigma_mat,
    zero_mat = para_list$zero_mat,
    quantile_mat = NULL,
    copula_list = copula_res$copula_list,
    n_cores = n_cores,
    family_use = family_use,
    nonnegative = nonnegative,
    nonzerovar = nonzerovar,
    input_data = input_data$dat,
    new_covariate = input_data$newCovariate,
    important_feature = copula_res$important_feature,
    parallelization = parallelization,
    BPPARAM = BPPARAM
  )
  message("New Data Generating End")

  scdesign3_res <- list(
    new_count = new_count,
    new_covariate = input_data$newCovariate,
    model_aic = copula_res$model_aic,
    model_bic = copula_res$model_bic,
    marginal_list = if (return_model)
      marginal_res
    else
      NULL,
    corr_list = if (return_model)
      copula_res$copula_list
    else
      NULL
  )
  return(scdesign3_res)
}

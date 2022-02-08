#' Simulate new data
#'
#' \code{simu_new} generates new simulated data based on fitted marginal and copula models.
#'
#' The function takes the new covariate (if use) from \code{\link{construct_data}},
#' parameter matricies from \code{\link{extract_para}} and multivariate Unifs from \code{\link{fit_copula}}.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param mean_mat A cell by feature matrix of the mean parameter.
#' @param sigma_mat A cell by feature matrix of the sigma parameter.
#' @param zero_mat A cell by feature matrix of the zero-inflation parameter.
#' @param quantile_mat A cell by feature matrix of the multivariate quantile.
#' @param n_cores An integer. The number of cores to use.
#' @param family_use A string of the marginal distribution.
#' Must be one of 'poisson', "binomial", 'nb', 'zip', 'zinb' or 'gaussian'.
#' @param nonnegative A logical variable. If TRUE, values < 0 will be converted to 0.Default is TRUE.
#' @param new_covariate A data.frame which contains covariates of targeted simulated data from  \code{\link{construct_data}}.
#'
#' @return A feature by cell matrix of the new simulated count (expression) matrix.
#'
#' @export simu_new

simu_new <- function(sce,
                     mean_mat,
                     sigma_mat,
                     zero_mat,
                     quantile_mat,
                     n_cores,
                     family_use,
                     nonnegative = TRUE,
                     new_covariate){
  ## New count
  mat <-  pbmcapply::pbmcmapply(function(x, y) {

    para_mat <- cbind(mean_mat[, x], sigma_mat[, x], quantile_mat[, x], zero_mat[, x])

    if (y == "binomial") {
      qfvec <- stats::qbinom(p = para_mat[, 3], prob = para_mat[, 1], size = 1)
    } else if (y == "poisson") {
      qfvec <- stats::qpois(p = para_mat[, 3], lambda = para_mat[, 1])
    } else if (y == "gaussian") {
      qfvec <-
        gamlss.dist::qNO(p = para_mat[, 3],
                         mu = para_mat[, 1],
                         sigma = abs(para_mat[, 2]))
    } else if (y == "nb") {
      qfvec <-
        gamlss.dist::qNBI(p = para_mat[, 3],
                          mu = para_mat[, 1],
                          sigma = para_mat[, 2])
    } else if (y == "zip") {
      qfvec <-
        gamlss.dist::qZIP(p = para_mat[, 3],
                          mu = para_mat[, 1],
                          sigma = para_mat[, 4])
    } else if (y == "zinb") {
      qfvec <-
        gamlss.dist::qZINBI(p = para_mat[, 3],
                            mu = para_mat[, 1],
                            sigma = para_mat[, 2],
                            nu = para_mat[, 4]
        )
    } else {
      stop("Distribution of gamlss must be one of gaussian, poisson, nb, zip or zinb!")
    }

    #message(paste0("Gene ", x , " End!"))

    r <- as.vector(qfvec)
    r
  }, x = seq_len(dim(sce)[1]), y = family_use, SIMPLIFY = TRUE, mc.cores = n_cores)

  new_count <- mat #simplify2array(mat)

  new_count <- as.matrix(t(new_count))

  if(nonnegative) new_count[new_count < 0] <- 0

  if(is.null(new_covariate)){
    colnames(new_count) <- colnames(sce)
    rownames(new_count) <- rownames(sce)
  }else{
    colnames(new_count) <- rownames(new_covariate)
    rownames(new_count) <- rownames(sce)
  }

  return(new_count)
}



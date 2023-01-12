#' Extract the parameters of each cell's distribution
#'
#' \code{extract_para} generates parameter matricies which determine each cell's distribution
#'
#' The function takes the new covariate (if use) from \code{\link{construct_data}} and
#' marginal models from \code{\link{fit_marginal}}.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param marginal_list A list of fitted regression models from \code{\link{fit_marginal}} for each gene in sce.
#' @param n_cores An integer. The number of cores to use.
#' @param family_use A string of the marginal distribution.
#' Must be one of 'poisson', 'nb', 'zip', 'zinb' or 'gaussian', which represent 'poisson distribution',
#' 'negative binomial distribution', 'zero-inflated poisson distribution', 'zero-inflated negative binomail distribution',
#' and 'gaussian distribution' respectively.
#' @param new_covariate A data.frame which contains covaraites of targeted simulated data from  \code{\link{construct_data}} and the
#' correlation group assignment for each cell in the column 'corr_group'.
#' @param parallelization A string indicating the specific parallelization function to use.
#' Must be one of 'mcmapply', 'bpmapply', or 'pbmcmapply', which corresponds to the parallelization function in the package
#' \code{parallel},\code{BiocParallel}, and \code{pbmcapply} respectively. The default value is 'mcmapply'.
#' @param BPPARAM A \code{MulticoreParam} object or NULL. When the parameter parallelization = 'mcmapply' or 'pbmcmapply',
#' this parameter must be NULL. When the parameter parallelization = 'bpmapply',  this parameter must be one of the
#' \code{MulticoreParam} object offered by the package 'BiocParallel. The default value is NULL.
#' @param data A dataframe which is used when fitting the gamlss model
#' @return A list with the components:
#' \describe{
#'   \item{\code{mean_mat}}{A cell by feature matrix of the mean parameter.}
#'   \item{\code{sigma_mat}}{A cell by feature matrix of the sigma parameter (for Gaussian, the variance; for NB, the dispersion.).}
#'   \item{\code{zero_mat}}{A cell by feature matrix of the zero-inflation parameter (only non-zero for ZIP and ZINB).}
#' }
#' @export extract_para

extract_para <-  function(sce,
                          marginal_list,
                          n_cores,
                          family_use,
                          new_covariate,
                          parallelization = "mcmapply",
                          BPPARAM = NULL,
                          data = NULL) {

  mat_function <-function(x, y) {
    fit <- marginal_list[[x]]

    if (methods::is(fit, "gamlss")) {
      mean_vec <-
        stats::predict(fit,
                       type = "response",
                       what = "mu",
                       newdata = new_covariate, data = data)
      if (y == "poisson" | y == "binomial") {
        theta_vec <- rep(NA, length(mean_vec))
      } else if (y == "gaussian") {
        theta_vec = stats::predict(fit,
                                   type = "response",
                                   what = "sigma",
                                   newdata = new_covariate, data = data) # this thete_vec is used for sigma_vec
      } else if (y == "nb") {
        theta_vec <-
          stats::predict(fit,
                         type = "response",
                         what = "sigma",
                         newdata = new_covariate, data = data)
        #theta_vec[theta_vec < 1e-3] <- 1e-3
      } else if (y == "zip") {
        theta_vec <- rep(NA, length(mean_vec))
        zero_vec <-
          stats::predict(fit,
                         type = "response",
                         what = "sigma",
                         newdata = new_covariate, data = data)
      } else if (y == "zinb") {
        theta_vec <-
          stats::predict(fit,
                         type = "response",
                         what = "sigma",
                         newdata = new_covariate, data = data)
        zero_vec <-
          stats::predict(fit,
                         type = "response",
                         what = "nu",
                         newdata = new_covariate, data = data)
      } else {
        stop("Distribution of gamlss must be one of gaussian, binomial, poisson, nb, zip or zinb!")
      }
    } else {
      ## Fit is mgcv::gam
      if (is.null(new_covariate)) {
        y <- stats::family(fit)$family[1]
        if (grepl("Negative Binomial", y)) {
          y <- "nb"
        }

        mean_vec <- stats::predict(fit, type = "response")
        if (y == "poisson" | y == "binomial") {
          theta_vec <- rep(NA, length(mean_vec))
        } else if (y == "gaussian") {
          theta_vec <- rep(sqrt(fit$sig2), length(mean_vec)) # this thete_vec is used for sigma_vec
        } else if (y == "nb") {
          theta <- fit$family$getTheta(TRUE)
          theta_vec <- 1/rep(theta, length(mean_vec))
        } else {
          stop("Distribution of gamlss must be one of gaussian, binomial, poisson, nb!")
        }
      } else{
        y <- stats::family(fit)$family[1]
        if (grepl("Negative Binomial", y)) {
          y <- "nb"
        }

        mean_vec <-
          stats::predict(fit, type = "response", newdata = new_covariate)
        if (y == "poisson" | y == "binomial") {
          theta_vec <- rep(NA, length(mean_vec))
        } else if (y == "gaussian") {
          theta_vec = stats::predict(fit,
                                     type = "response",
                                     what = "sigma",
                                     newdata = new_covariate) # this thete_vec is used for sigma_vec
        } else if (y == "nb") {
          theta <- fit$family$getTheta(TRUE)
          theta_vec <- 1/rep(theta, length(mean_vec))
        } else {
          stop("Distribution of gamlss must be one of gaussian, binomial, poisson, nb!")
        }
      }
    }

    #q_vec <- quantile_mat[, x]

    if (!exists("zero_vec")) {
      zero_vec <- 0
    }
    para_mat <- cbind(mean_vec, theta_vec, zero_vec)
    if (is.null(new_covariate)) {
      rownames(para_mat) <- colnames(sce)
    } else{
      rownames(para_mat) <- rownames(new_covariate)
    }
    para_mat
  }
  paraFunc <- parallel::mcmapply
  if(parallelization == "bpmapply"){
    paraFunc <- BiocParallel::bpmapply
  }
  if(parallelization == "pbmcmapply"){
    paraFunc <- pbmcapply::pbmcmapply
  }

  if(parallelization == "bpmapply"){
    BPPARAM$workers <- n_cores
    mat <- suppressMessages(paraFunc(mat_function, x = seq_len(dim(sce)[1]), y = family_use,BPPARAM = BPPARAM,SIMPLIFY = FALSE))
  }else{
    mat <- suppressMessages(paraFunc(mat_function, x = seq_len(dim(sce)[1]), y = family_use,SIMPLIFY = FALSE, mc.cores = n_cores))
  }
  mean_mat <- sapply(mat, function(x)
    x[, 1])
  sigma_mat <- sapply(mat, function(x)
    x[, 2])
  zero_mat <- sapply(mat, function(x)
    x[, 3])

  colnames(mean_mat) <-
    colnames(sigma_mat) <- colnames(zero_mat) <- rownames(sce)

  return(list(
    mean_mat = mean_mat,
    sigma_mat = sigma_mat,
    zero_mat = zero_mat
  ))
}

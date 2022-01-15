#' Simulate new data
#'
#' \code{simu_new} generates new simulated data based on fitted marginal and copula models.
#'
#' The function takes the new covariate (if use) from \code{\link{construct_data}},
#' marginal models from \code{\link{fit_marginal}} and multivariate Unifs from \code{\link{fit_copula}}.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param marginal_list A list of fitted regression models from \code{\link{fit_marginal}}.
#' @param copula_list A list of fitted copula models from \code{\link{fit_copula}}.
#' @param n_cores An integer. The number of cores to use.
#' @param family A string of the marginal distribution.
#' Must be one of 'poisson', 'nb', 'zip', 'zinb' or 'gaussian'.
#' @param new_covariate A data.frame which contains covaraites of targeted simulated data from  \code{\link{construct_data}}.
#'
#' @return A list with the components:
#' \describe{
#'   \item{\code{new_count}}{A matrix of the new simulated count (expression) matrix.}
#'   \item{\code{corr_list}}{A list of the fitted copula model. If using Gaussian copula, a list of correlation matrices; if vine, a list of vine objects.}
#'   \item{\code{model_aic}}{A vector of the marginal AIC and the copula AIC.}
#'   \item{\code{model_bic}}{A vector of the marginal BIC and the copula BIC.}
#' }
#'
#' @export simu_new
simu_new <- function(sce,
                     marginal_list,
                     copula_list,
                     n_cores,
                     family,
                     new_covariate){
  ## New count
  new.count <- cal_dist(sce = sce,
                       marginal_list = marginal_list,
                       quantile_mat = copula_list$new_mvu,
                       n_cores = n_cores, family = family,new_covariate = new_covariate)
  new.count <- as.matrix(t(new.count))
  if(is.null(new_covariate)){
    colnames(new.count) <- colnames(sce)
    rownames(new.count) <- rownames(sce)
  }else{
    colnames(new.count) <- rownames(new_covariate)
    rownames(new.count) <- rownames(sce)
  }

  return(list(new_count = new.count, model_aic = copula_list$model_aic, model_bic = copula_list$model_bic, corr_list = copula_list$corr_list))
}




## Calculate new data by quantiles
cal_dist <- function(sce,
                    marginal_list,
                    quantile_mat,
                    n_cores,
                    family,
                    new_covariate) {

  # n cell
  ncell <- dim(sce)[2]

  mat <- lapply(seq_len(dim(sce)[1]), function(x) {
    fit <- marginal_list[[x]]

    if (methods::is(fit, "gamlss")) {
      mean_vec <- stats::predict(fit, type = "response", what = "mu", newdata = new_covariate)
      if (family == "poisson") {
        theta_vec <- rep(NA, length(mean_vec))
      } else if (family == "gaussian") {
        theta_vec = stats::predict(fit, type = "response", what = "sigma", newdata = new_covariate) # this thete_vec is used for sigma_vec
      } else if (family == "nb") {
        theta_vec <- 1 / stats::predict(fit, type = "response", what = "sigma", newdata = new_covariate)
        theta_vec[theta_vec < 1e-3] <- 1e-3
      } else if (family == "zip") {
        theta_vec <- rep(NA, length(mean_vec))
        zero_vec <- stats::predict(fit, type = "response", what = "sigma", newdata = new_covariate)
      } else if (family == "zinb") {
        theta_vec <- stats::predict(fit, type = "response", what = "sigma", newdata = new_covariate)
        zero_vec <- stats::predict(fit, type = "response", what = "nu", newdata = new_covariate)
      } else {
        stop("Distribution of gamlss must be one of gaussian, poisson, nb, zip or zinb!")
      }
    } else { ## Fit is mgcv::gam
      if(is.null(new_covariate)){
        family <- stats::family(fit)$family[1]
        if(grepl("Negative Binomial", family)) {family <- "nb"}

        mean_vec <- stats::predict(fit, type = "response")
        if (family == "poisson") {
          theta_vec <- rep(NA, length(mean_vec))
        } else if (family == "gaussian") {
          theta_vec = stats::predict(fit, type = "response", what = "sigma") # this thete_vec is used for sigma_vec
        } else if (family == "nb") {
          theta <- fit$family$getTheta(TRUE)
          theta_vec <- rep(theta, length(mean_vec))
        } else {
          stop("Distribution of gamlss must be one of gaussian, poisson, nb!")
        }
      }else{
        family <- stats::family(fit)$family[1]
        if(grepl("Negative Binomial", family)) {family <- "nb"}

        mean_vec <- stats::predict(fit, type = "response", newdata = new_covariate)
        if (family == "poisson") {
          theta_vec <- rep(NA, length(mean_vec))
        } else if (family == "gaussian") {
          theta_vec = stats::predict(fit, type = "response", what = "sigma", newdata = new_covariate) # this thete_vec is used for sigma_vec
        } else if (family == "nb") {
          theta <- fit$family$getTheta(TRUE)
          theta_vec <- rep(theta, length(mean_vec))
        } else {
          stop("Distribution of gamlss must be one of gaussian, poisson, nb!")
        }
      }
    }

    q_vec <- quantile_mat[, x]

    if(!exists("zero_vec")) {zero_vec <- 0}
    para_mat <- cbind(mean_vec, theta_vec, q_vec, zero_vec)

    if (family == "poisson") {
      qfvec <- stats::qpois(p = para_mat[, 3], lambda = para_mat[, 1])
    } else if (family == "gaussian") {
      qfvec <-
        gamlss.dist::qNO(p = para_mat[, 3],
                         mu = para_mat[, 1],
                         sigma = abs(para_mat[, 2]))
    } else if (family == "nb") {
      qfvec <-
        gamlss.dist::qNBI(p = para_mat[, 3],
                          mu = para_mat[, 1],
                          sigma = 1 / para_mat[, 2])
    } else if (family == "zip") {
      qfvec <-
        gamlss.dist::qZIP(p = para_mat[, 3],
                          mu = para_mat[, 1],
                          sigma = para_mat[, 4])
    } else if (family == "zinb") {
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
  })

  new_mat <- simplify2array(mat)

  return(new_mat)
}

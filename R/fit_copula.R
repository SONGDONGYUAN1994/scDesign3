#' Fit the copula model
#'
#' \code{fit_copula} fits the copula model.
#'
#' This function takes the result from \code{\link{fit_marginal}} as the input and
#' and fit the copula model on the residuals.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param assay_use A string which indicates the assay you will use in the sce.
#' Default is 'counts'.
#' @param new_covariate A data.frame which contains covaraites of targeted simulated data from  \code{\link{construct_data}}.
#' @param marginal_list A list of fitted regression models from \code{\link{fit_marginal}}.
#' @param family A string of the marginal distribution. Must be one of 'poisson', 'nb', 'zip', 'zinb' or 'gaussian'.
#' @param copula A string of the copula choice. Must be one of 'gaussian' or 'vine'. Default is 'vine'.
#' @param DT A logic variable. If TRUE, perform the distributional transformation
#' to make the discrete data 'continuous'. This is useful for discrete distributions (e.g., Poisson, NB).
#' Default is TRUE.
#' @param pseudo_obs A logic variable. If TRUE, use the empirical quantiles instead of theoretical quantiles for fitting copula.
#' Default is FALSE.
#' @param epsilon A numeric variable for preventing the transformed quantiles to collapse to 0 or 1.
#' @param group A string or a string vector which indicates the groups for correlation structure.
#' @param ind A logic variable. If TRUE, no copula model will be fitted. Default is FASLE.
#' @param family_set A string or a string vector of the bivarate copula families. Default is c("gaussian", "indep").
#' @param n_cores An integer. The number of cores to use.
#'
#' @return A list with the components:
#' \describe{
#'   \item{\code{new_mvu}}{A matrix of the new multivariate uniform distribution from the copula.}
#'   \item{\code{cor_list}}{A list of the fitted copula model. If using Gaussian copula, a list of correlation matrices; if vine, a list of vine objects.}
#'   \item{\code{model_aic}}{A vector of the marginal AIC and the copula AIC.}
#' }
#'
#' @import mclust
#'
#' @export fit_copula

fit_copula <- function(sce,
                      assay_use,
                      new_covariate = NULL,
                      marginal_list,
                      family,
                      copula = 'vine',
                      DT = TRUE,
                      pseudo_obs = FALSE,
                      epsilon = 1e-6,
                      group,
                      ind = FALSE,
                      family_set = c("gaussian", "indep"),
                      n_cores) {
  # convert count matrix
  if (copula == "gaussian") {
    message("Convert Residuals to Multivariate Gaussian")
    newmat <- convert_n(
      sce = sce,
      assay_use = assay_use,
      marginal_list = marginal_list,
      DT = DT,
      pseudo_obs = pseudo_obs,
      family = family,
      epsilon = epsilon
    )
    message("Converting End")
  } else{
    message("Convert Residuals to Uniform")
    newmat <- convert_u(
      sce = sce,
      assay_use = assay_use,
      marginal_list = marginal_list,
      DT = DT,
      pseudo_obs = pseudo_obs,
      family = family,
      epsilon = epsilon
    )
    message("Converting End")
  }

  # identify groups
  ngene <- dim(sce)[1]
  ncell <- dim(sce)[2]

  if (is.null(rownames(newmat))) {
    rownames(newmat) <- colnames(sce)
  }

  group <- unlist(group)
  group_var <- unlist(group)

  if (group[1] == "none") {
    group <- rep(1, ncell)
  } else if (group[1] == "pseudotime" | length(group) > 1) {
    ## For continuous pseudotime, discretize it
    group <- SummarizedExperiment::colData(sce)[, group]
    mclust_mod <- mclust::Mclust(group, G = 1:5)

    group <- mclust_mod$classification

  } else {
    group <- SummarizedExperiment::colData(sce)[, group]
  }
  group_index <- unique(group)

  ## cluster for new covariate
  if(!is.null(new_covariate) & group_var[1] != "none"){
    if(group_var[1] == "pseudotime" | length(group_var) > 1){
      new_cov_group <- new_covariate[, group_var]
      mclust_mod <- mclust::Mclust(new_cov_group, G = 1:5) ## Here needs check. The group based on new covariate should correspond to group from ori covaraite.
      new_cov_group <- as.data.frame(mclust_mod$classification)
    }
  }
  # fit copula
  newmvn.list <- lapply(group_index, function(x, sce, newmat, ind, n_cores) {
    message(paste0("Copula group ", x, " starts"))
    curr_index <- colnames(sce)[which(group == x)]
    if(is.null(new_covariate)){
      curr_ncell <- length(curr_index)
      curr_ncell_idx <- curr_index
    }else{
      if(group_var[1] == "none"){
        curr_ncell <- dim(new_covariate)[1]
        curr_ncell_idx <- rownames(new_covariate)
      }else if(group_var[1] == "pseudotime" | length(group_var) > 1){
        curr_ncell <-  length(new_cov_group[which(new_cov_group==x),])
        curr_ncell_idx <- paste0("Cell", which(new_cov_group==x))
      }else{
        curr_ncell <-dim(as.data.frame(new_covariate[new_covariate[,group_var]==x,]))[1]
        curr_ncell_idx <- paste0("Cell", which(new_covariate[ ,1]==x))
      }
    }

    if (copula == "gaussian") {
      #message(paste0("Group ", group_index, " Start"))
      curr_mat <- newmat[curr_index, , drop = FALSE]
      #message("Cal MVN")
      cor.mat <- cal_cor(
        curr_mat,
        if.sparse = FALSE,
        lambda = 0.05,
        tol = 1e-8,
        ind = ind
      )
      #message("Sample MVN")
      new_mvu <- sampleMVN(n = curr_ncell,
                           Sigma = cor.mat)
      #message("MVN Sampling End")
      rownames(new_mvu) <- curr_ncell_idx

      #message("Cal AIC Start")
      model_aic <- cal_aic(norm.mat = newmat,
                          cor.mat = cor.mat,
                          ind = ind)
      #message("Cal AIC End")

    }else if(copula == "vine"){
      if(!ind) {
        curr_mat <- newmat[curr_index, , drop = FALSE]
        vine.fit <- rvinecopulib::vinecop(data = curr_mat,
                                          family_set = family_set,
                                          cores = n_cores)
        if(curr_ncell != 0){
          new_mvu <- rvinecopulib::rvinecop(curr_ncell,
                                            vine = vine.fit, cores = n_cores)
          rownames(new_mvu) <- curr_ncell_idx
        }else{
          new_mvu <- NULL
        }


        model_aic <- stats::AIC(vine.fit)

        cor.mat <- vine.fit
      }
      else {
        new_mvu <- matrix(data = stats::runif(curr_ncell*ngene), nrow = curr_ncell)
        rownames(new_mvu) <- curr_ncell_idx
        model_aic <- 0
        cor.mat <- NULL
      }
    }else{
      stop("Copula must be one of 'vine' or 'gaussian'")
    }
    return(list(
      new_mvu = new_mvu,
      model_aic = model_aic,
      cor.mat = cor.mat
    ))
  }, sce = sce, newmat = newmat, ind = ind, n_cores = n_cores)

  newmvn <- do.call(rbind, lapply(newmvn.list, function(x) x$new_mvu))

  copula.aic <- sum(sapply(newmvn.list, function(x) x$model_aic))

  marginal.aic <- sum(sapply(marginal_list, stats::AIC))

  model_aic <- c(marginal.aic, copula.aic)
  names(model_aic) <- c("aic.marginal", "aic.copula")

  cor_list <- lapply(newmvn.list, function(x) x$cor.mat)

  if(is.null(new_covariate)){
    new_mvu <- as.data.frame(newmvn[colnames(sce), ])
  }else{
    new_mvu <- as.data.frame(newmvn[rownames(new_covariate), ])
  }


  return(list(new_mvn = new_mvu, model_aic = model_aic, cor_list = cor_list))
}






## Calculate the correlation matrix. If use sparse cor estimation, package spcov will be used (it can be VERY SLOW).
cal_cor <- function(norm.mat,
                   if.sparse = FALSE,
                   lambda = 0.05,
                   tol = 1e-8,
                   ind = FALSE) {
  if (ind) {
    cor.mat <- diag(rep(1, dim(norm.mat)[2]))
    return(cor.mat)
  }
  else {
    cor.mat <- stats::cor(norm.mat, method = "pearson")
    s_d <- apply(norm.mat, 2, stats::sd)
    if (any(0 == s_d)) {
      cor.mat[is.na(cor.mat)] <- 0
    }
  }

  n <- dim(cor.mat)[1]

  if (if.sparse) {
    ifposd <- matrixcalc::is.positive.definite(cor.mat, tol = tol)
    if (!ifposd) {
      warning("Cor matrix is not positive defnite! Add tol to the diagnol.")
      #diag(cor.mat) <- diag(cor.mat) + tol
    }
    ## Call spcov
    lambda_matrix <- matrix(rep(1, n ^ 2), nrow = n) * lambda
    diag(lambda_matrix) <- 0
    scor <- spcov::spcov(
      cor.mat,
      cor.mat,
      lambda = lambda_matrix,
      step.size  = 100,
      trace = 1,
      n.inner.steps = 200,
      thr.inner = tol
    )
    cor.mat <- scor$Sigma
  }

  cor.mat
}

### Sample MVN based on cor
sampleMVN <- function(n,
                      Sigma) {
  mvnrv <-
    mvtnorm::rmvnorm(n, mean = rep(0, dim(Sigma)[1]), sigma = Sigma)
  mvnrvq <- apply(mvnrv, 2, stats::pnorm)

  return(mvnrvq)
}


## Convert marginal distributions to standard normals.
convert_n <- function(sce,
                     assay_use,
                     marginal_list,
                     DT = TRUE,
                     pseudo_obs = FALSE,
                     epsilon = 1e-6,
                     family) {
  ## Extract count matrix
  count_mat <- t(SummarizedExperiment::assay(sce, assay_use))

  # n cell
  ncell <- dim(count_mat)[1]

  mat <- sapply(seq_len(dim(sce)[1]), function(x) {
    fit <- marginal_list[[x]]

    if (methods::is(fit, "gamlss")) {
      mean_vec <- stats::predict(fit, type = "response", what = "mu")
      if (family == "poisson") {
        theta_vec <- rep(NA, length(mean_vec))
      } else if (family == "gaussian") {
        theta_vec <-
          stats::predict(fit, type = "response", what = "sigma") # called the theta_vec but actually used as sigma_vec for Gaussian
      } else if (family == "nb") {
        theta_vec <- 1 / stats::predict(fit, type = "response", what = "sigma")
        #theta_vec[theta_vec < 1e-3] <- 1e-3
      } else if (family == "zip") {
        theta_vec <- rep(NA, length(mean_vec))
        zero_vec <- stats::predict(fit, type = "response", what = "sigma")
      } else if (family == "zinb") {
        theta_vec <- stats::predict(fit, type = "response", what = "sigma")
        zero_vec <- stats::predict(fit, type = "response", what = "nu")
      } else {
        stop("Distribution of gamlss must be one of gaussian, poisson, nb, zip or zinb!")
      }
    } else {
      ## if input is from mgcv
      ## Check the family (since sometimes zip and zinb may degenerate into poisson or nb)

      family <- stats::family(fit)$family[1]
      if(grepl("Negative Binomial", family)) {family <- "nb"}

      mean_vec <- stats::predict(fit, type = "response")
      if (family == "poisson") {
        theta_vec <- rep(NA, length(mean_vec))
      } else if (family == "gaussian") {
        theta_vec <-
          stats::predict(fit, type = "response", what = "sigma") # called the theta_vec but actually used as sigma_vec for Gaussian
      } else if (family == "nb") {
        theta <- fit$family$getTheta(TRUE)
        theta_vec <- rep(theta, length(mean_vec))
      } else {
        stop("Distribution of mgcv must be one of gaussian, poisson or nb!")
      }
    }

    ## Mean for Each Cell


    Y <- count_mat[, x]


    ## Frame
    if(!exists("zero_vec")) {zero_vec <- 0}
    family_frame <- cbind(Y, mean_vec, theta_vec, zero_vec)

    if (family == "poisson") {
      pvec <- apply(family_frame, 1, function(x) {
        stats::ppois(x[1], lambda = x[2])
      })
    } else if (family == "gaussian") {
      pvec <- apply(family_frame, 1, function(x) {
        gamlss.dist::pNO(x[1], mu = x[2], sigma = abs(x[3]))
      })
    } else if (family == "nb") {
      pvec <- apply(family_frame, 1, function(x) {
        stats::pnbinom(x[1], mu = x[2], size = x[3])
      })
    } else if (family == "zip") {
      pvec <- apply(family_frame, 1, function(x) {
        gamlss.dist::pZIP(x[1], mu = x[2], sigma = abs(x[4]))
      })
    }
    else if (family == "zinb") {
      pvec <- apply(family_frame, 1, function(x) {
        gamlss.dist::pZINBI(x[1], mu = x[2], sigma = abs(x[3]), nu = x[4])
      })
    } else {
      stop("Distribution of gamlss must be one of gaussian, poisson, nb, zip or zinb!")
    }

    ## CHECK ABOUT THE FIRST PARAM!!!!!
    if (DT) {
      if (family == "poisson") {
        pvec2 <- apply(family_frame, 1, function(x) {
          stats::ppois(x[1] - 1, lambda = x[2]) * as.integer(x[1] > 0)
        })
      } else if (family == "gaussian") {
        ## Gaussian is continuous, thus do not need DT.
        message("Continuous gaussian does not need DT.")
        pvec2 <- pvec
        # pvec2 <- apply(family_frame, 1, function(x){
        # pNO(x[1]-1, mu = x[2], sigma = abs(x[3]))*as.integer(x[1] > 0)
        #})
      } else if (family == "nb") {
        pvec2 <- apply(family_frame, 1, function(x) {
          stats::pnbinom(x[1] - 1, mu = x[2], size = x[3]) * as.integer(x[1] > 0)
        })
      } else if (family == "zip") {
        pvec2 <- apply(family_frame, 1, function(x) {
          ifelse(x[1] > 0, gamlss.dist::pZIP(x[1] - 1, mu = x[2], sigma = abs(x[4])), 0)
        })
      } else if (family == "zinb") {
        pvec2 <- apply(family_frame, 1, function(x) {
          ifelse(x[1] > 0, gamlss.dist::pZINBI(x[1] - 1, mu = x[2], sigma = abs(x[3]), nu = x[4]), 0)
        })
      } else {
        stop("Distribution of gamlss must be one of gaussian, poisson, nb, zip or zinb!")
      }

      u1 <- pvec
      u2 <- pvec2

      v <- stats::runif(length(mean_vec))
      ## Random mapping
      r <- u1 * v + u2 * (1 - v)
    } else {
      r <- pvec
    }

    ## Avoid Inf
    idx_adjust <- which(1 - r < epsilon)
    r[idx_adjust] <- r[idx_adjust] - epsilon
    idx_adjust <- which(r < epsilon)
    r[idx_adjust] <- r[idx_adjust] + epsilon

    if (pseudo_obs) {
      r <- copula::pobs(r)
    }

    stats::qnorm(
      r,
      mean = 0,
      sd = 1,
      lower.tail = TRUE,
      log.p = FALSE
    )
  })

  colnames(mat) <- rownames(sce)
  rownames(mat) <- colnames(sce)

  ## Remove inf
  mat[is.infinite(mat)] <- NA

  for (i in 1:ncol(mat)) {
    mat[is.na(mat[, i]), i] <- mean(mat[, i], na.rm = TRUE)
  }

  return(mat)
}

## convert marginals to Unif[0, 1].
convert_u <- function(sce,
                     assay_use,
                     marginal_list,
                     DT = TRUE,
                     pseudo_obs = FALSE,
                     epsilon = 1e-6,
                     family) {
  ## Extract count matrix
  count_mat <- t(SummarizedExperiment::assay(sce, assay_use))

  # n cell
  ncell <- dim(count_mat)[1]

  mat <- sapply(seq_len(dim(sce)[1]), function(x) {
    fit <- marginal_list[[x]]

    if (methods::is(fit, "gamlss")) {
      mean_vec <- stats::predict(fit, type = "response", what = "mu")
      if (family == "poisson") {
        theta_vec <- rep(NA, length(mean_vec))
      } else if (family == "gaussian") {
        theta_vec <-
          stats::predict(fit, type = "response", what = "sigma") # called the theta_vec but actually used as sigma_vec for Gaussian
      } else if (family == "nb") {
        theta_vec <- 1 / stats::predict(fit, type = "response", what = "sigma")
        #theta_vec[theta_vec < 1e-3] <- 1e-3
      } else if (family == "zip") {
        theta_vec <- rep(NA, length(mean_vec))
        zero_vec <- stats::predict(fit, type = "response", what = "sigma")
      } else if (family == "zinb") {
        theta_vec <- stats::predict(fit, type = "response", what = "sigma")
        zero_vec <- stats::predict(fit, type = "response", what = "nu")
      } else {
        stop("Distribution of gamlss must be one of gaussian, poisson, nb, zip or zinb!")
      }
    } else {
      ## if input is from mgcv, update the family
      family <- stats::family(fit)$family[1]
      if(grepl("Negative Binomial", family)) {family <- "nb"}

      mean_vec <- stats::predict(fit, type = "response")
      if (family == "poisson") {
        theta_vec <- rep(NA, length(mean_vec))
      } else if (family == "gaussian") {
        theta_vec <-
          stats::predict(fit, type = "response", what = "sigma") # called the theta_vec but actually used as sigma_vec for Gaussian
      } else if (family == "nb") {
        theta <- fit$family$getTheta(TRUE)
        theta_vec <- rep(theta, length(mean_vec))
      } else {
        stop("Distribution of mgcv must be one of gaussian, poisson or nb!")
      }
    }

    ## Mean for Each Cell


    Y <- count_mat[, x]


    ## Frame
    if(!exists("zero_vec")) {zero_vec <- 0}
    family_frame <- cbind(Y, mean_vec, theta_vec, zero_vec)

    if (family == "poisson") {
      pvec <- apply(family_frame, 1, function(x) {
        stats::ppois(x[1], lambda = x[2])
      })
    } else if (family == "gaussian") {
      pvec <- apply(family_frame, 1, function(x) {
        gamlss.dist::pNO(x[1], mu = x[2], sigma = abs(x[3]))
      })
    } else if (family == "nb") {
      pvec <- apply(family_frame, 1, function(x) {
        stats::pnbinom(x[1], mu = x[2], size = x[3])
      })
    } else if (family == "zip") {
      pvec <- apply(family_frame, 1, function(x) {
        gamlss.dist::pZIP(x[1], mu = x[2], sigma = abs(x[4]))
      })
    }
    else if (family == "zinb") {
      pvec <- apply(family_frame, 1, function(x) {
        gamlss.dist::pZINBI(x[1], mu = x[2], sigma = abs(x[3]), nu = x[4])
      })
    } else {
      stop("Distribution of gamlss must be one of gaussian, poisson, nb, zip or zinb!")
    }

    ## CHECK ABOUT THE FIRST PARAM!!!!!
    if (DT) {
      if (family == "poisson") {
        pvec2 <- apply(family_frame, 1, function(x) {
          stats::ppois(x[1] - 1, lambda = x[2]) * as.integer(x[1] > 0)
        })
      } else if (family == "gaussian") {
        ## Gaussian is continuous, thus do not need DT.
        message("Continuous gaussian doesnot need DT.")
        pvec2 <- pvec
        # pvec2 <- apply(family_frame, 1, function(x){
        # pNO(x[1]-1, mu = x[2], sigma = abs(x[3]))*as.integer(x[1] > 0)
        #})
      } else if (family == "nb") {
        pvec2 <- apply(family_frame, 1, function(x) {
          stats::pnbinom(x[1] - 1, mu = x[2], size = x[3]) * as.integer(x[1] > 0)
        })
      } else if (family == "zip") {
        pvec2 <- apply(family_frame, 1, function(x) {
          ifelse(x[1] > 0, gamlss.dist::pZIP(x[1] - 1, mu = x[2], sigma = abs(x[4])), 0)
        })
      } else if (family == "zinb") {
        pvec2 <- apply(family_frame, 1, function(x) {
          ifelse(x[1] > 0, gamlss.dist::pZINBI(x[1] - 1, mu = x[2], sigma = abs(x[3]), nu = x[4]), 0)
        })
      } else {
        stop("Distribution of gamlss must be one of gaussian, poisson, nb, zip or zinb!")
      }

      u1 <- pvec
      u2 <- pvec2

      v <- stats::runif(length(mean_vec))
      ## Random mapping
      r <- u1 * v + u2 * (1 - v)
    } else {
      r <- pvec
    }

    ## Avoid Inf
    idx_adjust <- which(1 - r < epsilon)
    r[idx_adjust] <- r[idx_adjust] - epsilon
    idx_adjust <- which(r < epsilon)
    r[idx_adjust] <- r[idx_adjust] + epsilon

    if (pseudo_obs) {
      r <- copula::pobs(r)
    }

    r
  })

  colnames(mat) <- rownames(sce)
  rownames(mat) <- colnames(sce)

  ## Remove inf
  mat[is.infinite(mat)] <- NA

  ## Use mean to replace the missing value
  for (i in 1:ncol(mat)) {
    mat[is.na(mat[, i]), i] <- mean(mat[, i], na.rm = TRUE)
  }

  quantile_mat <- mat
  return(quantile_mat)
}


## Calcluate Gaussian copula aic
cal_aic <- function(norm.mat,
                   cor.mat,
                   ind) {

  if (ind) {
    copula.aic = 0
  } else {
    copula.nop <- (as.integer(sum(cor.mat != 0)) - dim(cor.mat)[1]) / 2

    copula.aic <-
      -2 * (sum(
        mvtnorm::dmvnorm(
          x = norm.mat,
          mean = rep(0, dim(cor.mat)[1]),
          sigma = cor.mat,
          log = TRUE
        )
      ) - sum(rowSums(stats::dnorm(norm.mat, log = TRUE)))) + 2 * copula.nop
  }

  copula.aic
}

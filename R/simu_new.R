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
#' @param copula_list A list of copulas for generating the multivariate quantile matrix. If provided, the \code{quantile_mat} must be NULL.
#' @param n_cores An integer. The number of cores to use.
#' @param family_use A string of the marginal distribution.
#' Must be one of 'poisson', "binomial", 'nb', 'zip', 'zinb' or 'gaussian'.
#' @param nonnegative A logical variable. If TRUE, values < 0 will be converted to 0. Default is TRUE (since the expression matrix is nonnegative).
#' @param nonzerovar A logical variable. If TRUE, for any gene with zero variance, a cell will be replaced with 1. This is designed for avoiding potential errors, for example, PCA.
#' @param input_data A input count matrix.
#' @param new_covariate A data.frame which contains covariates of targeted simulated data from  \code{\link{construct_data}}.
#'
#' @return A feature by cell matrix of the new simulated count (expression) matrix.
#'
#' @export simu_new

simu_new <- function(sce,
                     mean_mat,
                     sigma_mat,
                     zero_mat,
                     quantile_mat = NULL,
                     copula_list,
                     n_cores,
                     family_use,
                     nonnegative = TRUE,
                     nonzerovar = TRUE,
                     input_data,
                     new_covariate){
  if(!is.null(quantile_mat) & !is.null(copula_list)) {
    stop("You can only provide either the quantile_mat or the copula_list!")
  }

  if(!is.null(quantile_mat)) {
    message("Multivariate quantile matrix is provided")
  } else {
    message("Use Copula to sample a multivariate quantile matrix")

    group_index <- unique(input_data$corr_group)
    corr_group <- as.data.frame(input_data$corr_group)
    colnames(corr_group) <- "corr_group"
    ngene <- dim(sce)[1]
    if (is.null(new_covariate)) {
      new_corr_group <- NULL
    } else{
      new_corr_group <- as.data.frame(new_covariate$corr_group)
      colnames(new_corr_group) <- "corr_group"
    }
    ind <- group_index[1] == "ind"
    newmvn.list <-
      lapply(group_index, function(x,
                                   sce,
                                   corr_group,
                                   new_corr_group,
                                   ind,
                                   n_cores,
                                   copula_list) {
        message(paste0("Sample Copula group ", x, " starts"))
        curr_index <- which(corr_group[, 1] == x)
        if (is.null(new_covariate)) {
          curr_ncell <- length(curr_index)
          curr_ncell_idx <- curr_index
        } else{
          curr_ncell <- length(which(new_corr_group[, 1] == x))
          curr_ncell_idx <-
            paste0("Cell", which(new_corr_group[, 1] == x))
        }
        cor.mat <- copula_list[[x]]

        if(curr_ncell == 0) {
          new_mvu <- NULL
        } else {
          if (class(cor.mat)[1] == "matrix") {
            #message(paste0("Group ", group_index, " Start"))

            #message("Sample MVN")
            new_mvu <- sampleMVN(n = curr_ncell,
                                 Sigma = cor.mat)
            #message("MVN Sampling End")
            rownames(new_mvu) <- curr_ncell_idx
          } else if (class(cor.mat)[1] == "vinecop") {

            #message("Sampling Vine Copula Starts")
            new_mvu <- rvinecopulib::rvinecop(
              curr_ncell,
              vine = cor.mat,
              cores = n_cores,
              qrng = TRUE
            )
            #message("Sampling Vine Copula Ends")
            rownames(new_mvu) <- curr_ncell_idx
          } else if (ind) {
            "Use independent copula (random Unif)."
            new_mvu <-
              matrix(data = stats::runif(curr_ncell * ngene),
                     nrow = curr_ncell)
            rownames(new_mvu) <- curr_ncell_idx
          } else{
            stop("Copula must be one from 'vine' or 'gaussian', or assume gene-gene is independent")
          }
        }
        return(
          list(
            new_mvu = new_mvu
          )
        )
      }, sce = sce, ind = ind, n_cores = n_cores, corr_group = corr_group, new_corr_group = new_corr_group, copula_list = copula_list)

    newmvn <-
      do.call(rbind, lapply(newmvn.list, function(x)
        x$new_mvu))

    quantile_mat <- as.data.frame(newmvn)
  }

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
                          sigma = max(para_mat[, 4], 2.2e-16)) ## Avoid zero zero-inflated prob
    } else if (y == "zinb") {
      qfvec <-
        gamlss.dist::qZINBI(p = para_mat[, 3],
                            mu = para_mat[, 1],
                            sigma = para_mat[, 2],
                            nu = max(para_mat[, 4], 2.2e-16)
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

  if(nonzerovar) {
    row_vars <- Rfast::rowVars(new_count)
    if(sum(row_vars == 0) > 0) {
      message("Some genes have zero variance. Replace a random one with 1.")
      row_vars_index <- which(row_vars == 0)
      col_index <- seq_len(dim(new_count)[2])
      for(i in row_vars_index) {
        new_count[i, sample(col_index, 1)] <- 1
      }
    }
  }

  return(new_count)
}



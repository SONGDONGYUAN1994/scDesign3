#' Extract the parameters of each cell's distribution
#'
#' \code{extract_para} generates parameter matricies which determine each cell's distribution
#'
#' The function takes the new covariate (if use) from \code{\link{construct_data}} and
#' marginal models from \code{\link{fit_marginal}}.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param assay_use A string which indicates the assay you will use in the sce. Default is 'counts'.
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
#' @examples
#'   data(example_sce)
#'   my_data <- construct_data(
#'   sce = example_sce,
#'   assay_use = "counts",
#'   celltype = "cell_type",
#'   pseudotime = "pseudotime",
#'   spatial = NULL,
#'   other_covariates = NULL,
#'   corr_by = "1"
#'   )
#'   my_marginal <- fit_marginal(
#'   data = my_data,
#'   mu_formula = "s(pseudotime, bs = 'cr', k = 10)",
#'   sigma_formula = "1",
#'   family_use = "nb",
#'   n_cores = 1,
#'   usebam = FALSE
#'   )
#'   my_copula <- fit_copula(
#'   sce = example_sce,
#'   assay_use = "counts",
#'   marginal_list = my_marginal,
#'   family_use = c(rep("nb", 5), rep("zip", 5)),
#'   copula = "vine",
#'   n_cores = 1,
#'   new_covariate = NULL,
#'   input_data = my_data$dat
#'   )
#'   my_para <- extract_para(
#'     sce = example_sce,
#'     marginal_list = my_marginal,
#'     n_cores = 1,
#'     family_use = c(rep("nb", 5), rep("zip", 5)),
#'     new_covariate = NULL,
#'     data = my_data$dat
#'   )
#'
#' @export extract_para

extract_para <-  function(sce,
                          assay_use = "counts",
                          marginal_list,
                          n_cores,
                          family_use,
                          new_covariate,
                          parallelization = "mcmapply",
                          BPPARAM = NULL,
                          data) {
  removed_cell_list <- lapply(marginal_list, function(x){x$removed_cell})
  marginal_list <- lapply(marginal_list, function(x){x$fit})

  # find gene whose marginal is fitted
  qc_gene_idx <- which(!is.na(marginal_list))

  mat_function <-function(x, y) {
    fit <- marginal_list[[x]]
    removed_cell <- removed_cell_list[[x]]
    count_mat <-
      t(as.matrix(SummarizedExperiment::assay(sce, assay_use)))
    data$gene <- count_mat[,x]
    # if(!"gamlss" %in% class(fit)){
    #   modelframe <- model.frame(fit)
    # }else{
    #   modelframe <- fit$mu.x
    # }
    if(is.null(new_covariate)){
      total_cells <- dim(sce)[2]
      cell_names <- colnames(sce)
    }else{
      total_cells <- dim(new_covariate)[1]
      cell_names <- rownames(new_covariate)
    }

    if(length(removed_cell) > 0 && !any(is.na(removed_cell))){
      if(is.null(new_covariate)){
        data <- data[-removed_cell,]
      }else{
        if (methods::is(fit, "gamlss")){
          all_covariates <- all.vars(fit$mu.formula)[-1]
        }else{
          all_covariates <- all.vars(fit$formula)[-1]
        }
        remove_idx <- lapply(all_covariates, function(x){
          curr_x <- tapply(data$gene, data[,x], sum)
          zero_group <- which(curr_x==0)
          if(length(zero_group) == 0){
            return(NA)
          }else{
            type <- names(curr_x)[zero_group]
            return(which(new_covariate[,x] %in% type))
          }

        })
        remove_cell_idx <- unlist(remove_idx)
        remove_cell_idx <- unique(stats::na.omit(remove_cell_idx))
        if(length(remove_cell_idx) > 0){
          new_covariate <- new_covariate[-remove_cell_idx,]
        }
      }
    }

    if (methods::is(fit, "gamlss")) {
      mean_vec <-
        stats::predict(fit,
                       type = "response",
                       what = "mu",
                       newdata = new_covariate, data = data)
      if (y == "poisson" | y == "binomial") {
        theta_vec <- rep(NA, total_cells)
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
        theta_vec <- rep(NA,  total_cells)
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
          theta_vec <- rep(NA,  total_cells)
        } else if (y == "gaussian") {
          theta_vec <- rep(sqrt(fit$sig2), total_cells) # this thete_vec is used for sigma_vec
        } else if (y == "nb") {
          theta <- fit$family$getTheta(TRUE)
          theta_vec <- 1/rep(theta, total_cells)
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
          theta_vec <- rep(NA, total_cells)
        } else if (y == "gaussian") {
          theta_vec = stats::predict(fit,
                                     type = "response",
                                     what = "sigma",
                                     newdata = new_covariate) # this thete_vec is used for sigma_vec
        } else if (y == "nb") {
          theta <- fit$family$getTheta(TRUE)
          theta_vec <- 1/rep(theta, total_cells)
        } else {
          stop("Distribution of gam must be one of gaussian, binomial, poisson, nb!")
        }
      }


    }

    if (!exists("zero_vec")) {
      zero_vec <- rep(0, length(mean_vec))
      names(zero_vec) <- names(mean_vec)
    }

    if(length(mean_vec) < total_cells){
      full_means <- rep(NA, total_cells)
      names(full_means) <- cell_names
      full_means[names(mean_vec)] <- mean_vec
      full_theta <- rep(NA, total_cells)
      names(full_theta) <- cell_names
      full_zero <- rep(NA, total_cells)
      names(full_zero) <- cell_names
      if(is.null(names(theta_vec))){
        if(length(theta_vec) == length(mean_vec)){
          names(theta_vec) <- names(mean_vec)  ## for gamlss case
        }else{
          names(theta_vec) <- cell_names ## for gam case
        }

      }
      full_theta[names(theta_vec)] <- theta_vec
      full_zero[names(zero_vec)] <- zero_vec
      mean_vec <- full_means
      theta_vec <- full_theta
      zero_vec <- full_zero
    }

    #q_vec <- quantile_mat[, x]


    para_mat <- cbind(mean_vec, theta_vec, zero_vec)
    rownames(para_mat) <- cell_names
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
    mat <- suppressMessages(paraFunc(mat_function, x = seq_len(dim(sce)[1])[qc_gene_idx], y = family_use,BPPARAM = BPPARAM,SIMPLIFY = FALSE))
  }else{
    mat <- suppressMessages(paraFunc(mat_function, x = seq_len(dim(sce)[1])[qc_gene_idx], y = family_use,SIMPLIFY = FALSE
                                   ,mc.cores = n_cores
                                   ))
  }
  mean_mat <- sapply(mat, function(x)
    x[, 1])
  sigma_mat <- sapply(mat, function(x)
    x[, 2])
  zero_mat <- sapply(mat, function(x)
    x[, 3])

  if(length(qc_gene_idx) > 0){
    colnames(mean_mat) <-
      colnames(sigma_mat) <- colnames(zero_mat) <- rownames(sce)[qc_gene_idx]
    na_mat <- matrix(NA, nrow = dim(mean_mat)[1], ncol = dim(sce)[1] - length(qc_gene_idx))
    rownames(na_mat) <- rownames(mean_mat)
    colnames(na_mat) <- rownames(sce)[-qc_gene_idx]
    mean_mat <- cbind(mean_mat,na_mat)
    sigma_mat <- cbind(sigma_mat,na_mat)
    zero_mat <- cbind(zero_mat,na_mat)
    mean_mat <- mean_mat[,rownames(sce)]
    sigma_mat <- sigma_mat[,rownames(sce)]
    zero_mat <- zero_mat[,rownames(sce)]
  }else{
    colnames(mean_mat) <-
      colnames(sigma_mat) <- colnames(zero_mat) <- rownames(sce)
  }


  return(list(
    mean_mat = mean_mat,
    sigma_mat = sigma_mat,
    zero_mat = zero_mat
  ))
}

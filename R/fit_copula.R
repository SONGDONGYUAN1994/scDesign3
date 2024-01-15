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
#' @param input_data The input data, which is one of the output from \code{\link{construct_data}}.
#' @param empirical_quantile Please only use it if you clearly know what will happen! A logic variable. If TRUE, DO NOT fit the copula and use the EMPIRICAL CDF values of the original data; it will make the simulated data fixed (no randomness). Default is FALSE. Only works if ncell is the same as your original data.
#' @param marginal_list A list of fitted regression models from \code{\link{fit_marginal}}.
#' @param family_use A string or a vector of strings of the marginal distribution. Must be one of 'poisson', 'nb', 'zip', 'zinb' or 'gaussian'.
#' @param copula A string of the copula choice. Must be one of 'gaussian' or 'vine'. Default is 'gaussian'. Note that vine copula may have better modeling of high-dimensions, but can be very slow when features are >1000.
#' @param DT A logic variable. If TRUE, perform the distributional transformation
#' to make the discrete data 'continuous'. This is useful for discrete distributions (e.g., Poisson, NB).
#' Default is TRUE. Note that for continuous data (e.g., Gaussian), DT does not make sense and should be set as FALSE.
#' @param pseudo_obs A logic variable. If TRUE, use the empirical quantiles instead of theoretical quantiles for fitting copula.
#' Default is FALSE.
#' @param epsilon A numeric variable for preventing the transformed quantiles to collapse to 0 or 1.
#' @param family_set A string or a string vector of the bivariate copula families. Default is c("gaussian", "indep").
#' @param important_feature A string or vector which indicates whether a gene will be used in correlation estimation or not. If this is a string, then
#' this string must be either "all" (using all genes) or "auto", which indicates that the genes will be automatically selected based on the proportion of zero expression across cells
#' for each gene. Gene with zero proportion greater than 0.8 will be excluded form gene-gene correlation estimation. If this is a vector, then this should
#' be a logical vector with length equal to the number of genes in \code{sce}. \code{TRUE} in the logical vector means the corresponding gene will be included in
#' gene-gene correlation estimation and \code{FALSE} in the logical vector means the corresponding gene will be excluded from the gene-gene correlation estimation.
#' The default value for is "all".
#' @param n_cores An integer. The number of cores to use.
#' @param parallelization A string indicating the specific parallelization function to use.
#' Must be one of 'mcmapply', 'bpmapply', or 'pbmcmapply', which corresponds to the parallelization function in the package
#' \code{parallel},\code{BiocParallel}, and \code{pbmcapply} respectively. The default value is 'mcmapply'.
#' @param BPPARAM A \code{MulticoreParam} object or NULL. When the parameter parallelization = 'mcmapply' or 'pbmcmapply',
#' this parameter must be NULL. When the parameter parallelization = 'bpmapply',  this parameter must be one of the
#' \code{MulticoreParam} object offered by the package 'BiocParallel. The default value is NULL.
#'
#' @return A list with the components:
#' \describe{
#'   \item{\code{new_mvu}}{A matrix of the new multivariate uniform distribution from the copula.}
#'   \item{\code{copula_list}}{A list of the fitted copula model. If using Gaussian copula, a list of correlation matrices; if vine, a list of vine objects.}
#'   \item{\code{model_aic}}{A vector of the marginal AIC and the copula AIC.}
#'   \item{\code{model_bic}}{A vector of the marginal BIC and the copula BIC.}
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
#'   input_data = my_data$dat
#'   )
#'   
#' @import mclust
#' @import gamlss
#'
#' @export fit_copula

fit_copula <- function(sce,
                       assay_use,
                       input_data,
                       empirical_quantile = FALSE,
                       marginal_list,
                       family_use,
                       copula = 'gaussian',
                       DT = TRUE,
                       pseudo_obs = FALSE,
                       epsilon = 1e-6,
                       family_set = c("gaussian", "indep"),
                       important_feature = "all",
                       n_cores,
                       parallelization = "mcmapply",
                       BPPARAM = NULL) {

  if(empirical_quantile == TRUE) {
      message("Use the empirical quantile matrices from the original data; do not fit copula. This will make the result FIXED.")
  }
  
  if(important_feature == "all") {
    important_feature <- rep(TRUE, dim(sce)[1])
  }
  
  marginals <- lapply(marginal_list, function(x){x$fit})
  # find gene whose marginal is fitted
  qc_gene_idx <- which(!is.na(marginals))
  if(length(family_use) != 1){
    family_use <- family_use[qc_gene_idx]
  }
  group_index <- unique(input_data$corr_group)
  ind <- group_index[1] == "ind"
  if(ind) {
    copula.aic <- 0
    copula.bic <- 0
    marginal.aic <- sum(sapply(marginals[qc_gene_idx], stats::AIC))
    marginal.bic <- sum(sapply(marginals[qc_gene_idx], stats::BIC))
    model_aic <- c(marginal.aic, copula.aic, marginal.aic + copula.aic)
    names(model_aic) <- c("aic.marginal", "aic.copula", "aic.total")
    model_bic <- c(marginal.bic, copula.bic, marginal.bic + copula.bic)
    names(model_bic) <- c("bic.marginal", "bic.copula", "bic.total")
    copula_list <- NULL
    return(
      list(
        #new_mvu = new_mvu,
        model_aic = model_aic,
        model_bic = model_bic,
        copula_list = copula_list,
        important_feature = important_feature
      )
    )
  }

  if(empirical_quantile == TRUE) {
    ###
    if(is.vector(important_feature) & methods::is(important_feature,"logical")){
      if(length(important_feature) != dim(sce)[1]){
        stop("The important_feature should either be 'auto' or a logical vector with the length equals to the number of genes in the input data")
      }
    }else{
      if(important_feature=="auto"){
        gene_zero_prop <- apply(as.matrix(SummarizedExperiment::assay(sce, assay_use)), 1, function(y){
          sum(y < 1e-5) / dim(sce)[2]
        })
        important_feature = gene_zero_prop < 0.8 ## default zero proportion in scDesign2
        names(important_feature) <- rownames(sce)
      }else{
        stop("The important_feature should either be 'auto' or a logical vector with the length equals to the number of genes in the input data")
      }
    }
    
    important_feature <- important_feature[qc_gene_idx]
    corr_group <- as.data.frame(input_data$corr_group)
    colnames(corr_group) <- "corr_group"

    
    newmvq.list <- lapply(group_index, function(x,
                                                sce,
                                                corr_group) {
      message(paste0("Empirical quantile group ", x, " starts"))
      sce <- sce[important_feature, ]
      curr_index <- which(corr_group[, 1] == x)
      
      newmat <- SummarizedExperiment::assay(sce[, curr_index], assay_use)
      newmat <- t(as.matrix(newmat))
      
      newmvq <- rvinecopulib::pseudo_obs(newmat)
      newmvq
    }, sce = sce,
    corr_group = corr_group)
      
    newmvq <- do.call("rbind", newmvq.list)
    newmvq <- newmvq[colnames(sce),] 
    
    return(list(model_aic = 0,
           model_bic = 0,
           quantile_mat = newmvq,
           important_feature = important_feature))
    ###
  } else {
    if (copula == "gaussian") {
      message("Convert Residuals to Multivariate Gaussian")
      newmat <- convert_n(
        sce = sce[qc_gene_idx,],
        assay_use = assay_use,
        marginal_list = marginal_list[qc_gene_idx],
        data = input_data,
        DT = DT,
        pseudo_obs = pseudo_obs,
        n_cores = n_cores,
        family_use = family_use,
        epsilon = epsilon,
        parallelization = parallelization,
        BPPARAM = BPPARAM
      )
      message("Converting End")
    } else{
      message("Convert Residuals to Uniform")
      newmat <- convert_u(
        sce = sce[qc_gene_idx,],
        assay_use = assay_use,
        marginal_list = marginal_list[qc_gene_idx],
        data = input_data,
        DT = DT,
        pseudo_obs = pseudo_obs,
        family_use = family_use,
        n_cores = n_cores,
        epsilon = epsilon,
        parallelization = parallelization,
        BPPARAM = BPPARAM
      )
      message("Converting End")
    }
    
    ## select important genes
    if(is.vector(important_feature) & methods::is(important_feature,"logical")){
      if(length(important_feature) != dim(sce)[1]){
        stop("The important_feature should either be 'auto' or a logical vector with the length equals to the number of genes in the input data")
      }
    }else{
      if(important_feature=="auto"){
        gene_zero_prop <- apply(as.matrix(SummarizedExperiment::assay(sce, assay_use)), 1, function(y){
          sum(y < 1e-5) / dim(sce)[2]
        })
        important_feature = gene_zero_prop < 0.8 ## default zero proportion in scDesign2
        names(important_feature) <- rownames(sce)
      }else{
        stop("The important_feature should either be 'auto' or a logical vector with the length equals to the number of genes in the input data")
      }
    }
    
    important_feature <- important_feature[qc_gene_idx]
    
    
    corr_group <- as.data.frame(input_data$corr_group)
    colnames(corr_group) <- "corr_group"
   
    
    newmvn.list <-
      lapply(group_index, function(x,
                                   sce,
                                   newmat,
                                   corr_group,
                                   ind,
                                   n_cores,
                                   important_feature) {
        message(paste0("Copula group ", x, " starts"))
        curr_index <- which(corr_group[, 1] == x)
       
        if (copula == "gaussian") {
          #message(paste0("Group ", group_index, " Start"))
          curr_mat <- newmat[curr_index, , drop = FALSE]
          #message("Cal MVN")
          cor.mat <- cal_cor(
            curr_mat,
            important_feature = important_feature,
            if.sparse = FALSE,
            lambda = 0.05,
            tol = 1e-8,
            ind = ind
          )
          
          #message("Sample MVN")
          #new_mvu <- sampleMVN(n = curr_ncell,
          #                     Sigma = cor.mat)
          #message("MVN Sampling End")
          #rownames(new_mvu) <- curr_ncell_idx
          
          #message("Cal AIC/BIC Start")
          model_aic <- cal_aic(norm.mat = newmat,
                               cor.mat = cor.mat,
                               ind = ind)
          model_bic <- cal_bic(norm.mat = newmat,
                               cor.mat = cor.mat,
                               ind = ind)
          #message("Cal AIC/BIC End")
          
        } else if (copula == "vine") {
          message("Vine Copula Estimation Starts")
          start <- Sys.time()
          curr_mat <- newmat[curr_index, , drop = FALSE]
          curr_mat <- curr_mat[,which(important_feature)]
          vine.fit <- rvinecopulib::vinecop(
            data = curr_mat,
            family_set = family_set,
            show_trace = FALSE,
            par_method = "mle",
            cores = n_cores
          )
          end <- Sys.time()
          print(end - start)
          message("Vine Copula Estimation Ends")
          # if (curr_ncell != 0) {
          #   message("Sampling Vine Copula Starts")
          #   new_mvu <- rvinecopulib::rvinecop(
          #     curr_ncell,
          #     vine = vine.fit,
          #     cores = n_cores,
          #     qrng = TRUE
          #   )
          #   message("Sampling Vine Copula Ends")
          #   rownames(new_mvu) <- curr_ncell_idx
          # } else{
          #   new_mvu <- NULL
          # }
          model_aic <- stats::AIC(vine.fit)
          model_bic <- stats::BIC(vine.fit)
          cor.mat <- vine.fit
        } else{
          stop("Copula must be one of 'vine' or 'gaussian'")
        }
        return(
          list(
            #new_mvu = new_mvu,
            model_aic = model_aic,
            model_bic = model_bic,
            cor.mat = cor.mat
          )
        )
      }, sce = sce, 
      newmat = newmat, 
      ind = ind, 
      n_cores = n_cores, 
      corr_group = corr_group,
      important_feature = important_feature)
    
    #newmvn <-
    #  do.call(rbind, lapply(newmvn.list, function(x)
    #    x$new_mvu))
    
    copula.aic <- sum(sapply(newmvn.list, function(x)
      x$model_aic))
    marginal.aic <- sum(sapply(marginals[qc_gene_idx], stats::AIC))
    
    copula.bic <- sum(sapply(newmvn.list, function(x)
      x$model_bic))
    marginal.bic <- sum(sapply(marginals[qc_gene_idx], stats::BIC))
    
    model_aic <- c(marginal.aic, copula.aic, marginal.aic + copula.aic)
    names(model_aic) <- c("aic.marginal", "aic.copula", "aic.total")
    
    model_bic <- c(marginal.bic, copula.bic, marginal.bic + copula.bic)
    names(model_bic) <- c("bic.marginal", "bic.copula", "bic.total")
    
    copula_list <- lapply(newmvn.list, function(x)
      x$cor.mat)
    names(copula_list) <- group_index
    
    
    #new_mvu <- as.data.frame(newmvn)
    
    return(
      list(
        #new_mvu = new_mvu,
        model_aic = model_aic,
        model_bic = model_bic,
        copula_list = copula_list,
        important_feature = important_feature
      )
    )
  }
}


## Calculate the correlation matrix. If use sparse cor estimation, package spcov will be used (it can be VERY SLOW).
cal_cor <- function(norm.mat,
                    important_feature,
                    if.sparse = FALSE,
                    lambda = 0.05,
                    tol = 1e-8,
                    ind = FALSE) {
  if (ind) {
    cor.mat <- diag(rep(1, dim(norm.mat)[2]))
    return(cor.mat)
  }
  else {
    cor.mat <- diag(rep(1, dim(norm.mat)[2]))
    rownames(cor.mat) <- colnames(norm.mat)
    colnames(cor.mat) <- colnames(norm.mat)
    important.mat <- norm.mat[,which(important_feature)]
    important_cor.mat <- corrlation(important.mat)
    #s_d <- apply(norm.mat, 2, stats::sd)
    s_d <- matrixStats::colSds(important.mat,na.rm = TRUE)
    if (any(0 == s_d)) {
      important_cor.mat[is.na(important_cor.mat)] <- 0
    }
    cor.mat[rownames(important_cor.mat), colnames(important_cor.mat)] <- important_cor.mat
  }

  n <- dim(cor.mat)[1]

  ### We currently gave up this because the sparse corr calclulation is too slow.
  # if (if.sparse) {
  #   ifposd <- matrixcalc::is.positive.definite(cor.mat, tol = tol)
  #   if (!ifposd) {
  #     warning("Cor matrix is not positive defnite! Add tol to the diagnol.")
  #     #diag(cor.mat) <- diag(cor.mat) + tol
  #   }
  #   ## Call spcov
  #   lambda_matrix <- matrix(rep(1, n ^ 2), nrow = n) * lambda
  #   diag(lambda_matrix) <- 0
  #   scor <- spcov::spcov(
  #     cor.mat,
  #     cor.mat,
  #     lambda = lambda_matrix,
  #     step.size  = 100,
  #     trace = 1,
  #     n.inner.steps = 200,
  #     thr.inner = tol
  #   )
  #   cor.mat <- scor$Sigma
  # }

  cor.mat
}

## Convert marginal distributions to standard normals.
convert_n <- function(sce,
                      assay_use,
                      marginal_list,
                      data,
                      DT = TRUE,
                      pseudo_obs = FALSE,
                      epsilon = 1e-6,
                      family_use,
                      n_cores,
                      parallelization,
                      BPPARAM) {
  ## Extract count matrix
  count_mat <-
      t(as.matrix(SummarizedExperiment::assay(sce, assay_use)))
  removed_cell_list <- lapply(marginal_list, function(x){x$removed_cell})
  marginal_list <- lapply(marginal_list, function(x){x$fit})
  # n cell
  ncell <- dim(count_mat)[1]


  mat_function <- function(x, y) {
    fit <- marginal_list[[x]]
    removed_cell <- removed_cell_list[[x]]
    if(length(removed_cell) > 0 && !any(is.na(removed_cell))){
      data<- data[-removed_cell,]
    }
    if (methods::is(fit, "gamlss")) {
      mean_vec <- stats::predict(fit, type = "response", what = "mu", data = data)
      if (y == "poisson" | y == "binomial") {
        theta_vec <- rep(NA, length(mean_vec))
      } else if (y == "gaussian") {
        theta_vec <-
          stats::predict(fit, type = "response", what = "sigma", data = data) # called the theta_vec but actually used as sigma_vec for Gaussian
      } else if (y == "nb") {
        theta_vec <-
          1 / stats::predict(fit, type = "response", what = "sigma", data = data)
        #theta_vec[theta_vec < 1e-3] <- 1e-3
      } else if (y == "zip") {
        theta_vec <- rep(NA, length(mean_vec))
        zero_vec <-
          stats::predict(fit, type = "response", what = "sigma", data = data)
      } else if (y == "zinb") {
        theta_vec <- stats::predict(fit, type = "response", what = "sigma", data = data)
        zero_vec <-
          stats::predict(fit, type = "response", what = "nu", data = data)
      } else {
        stop("Distribution of gamlss must be one of gaussian, poisson, nb, zip or zinb!")
      }
    } else {
      ## if input is from mgcv
      ## Check the family (since sometimes zip and zinb may degenerate into poisson or nb)

      y <- stats::family(fit)$family[1]
      if (grepl("Negative Binomial", y)) {
        y <- "nb"
      }

      mean_vec <- stats::predict(fit, type = "response")
      if (y == "poisson" | y == "binomial") {
        theta_vec <- rep(NA, length(mean_vec))
      } else if (y == "gaussian") {
        theta_vec <- rep(sqrt(fit$sig2), length(mean_vec)) # called the theta_vec but actually used as sigma_vec for Gaussian
      } else if (y == "nb") {
        theta <- fit$family$getTheta(TRUE)
        theta_vec <- rep(theta, length(mean_vec))
      } else {
        stop("Distribution of mgcv must be one of gaussian, poisson or nb!")
      }
    }

    ## Mean for Each Cell


    Y <- count_mat[names(mean_vec), x]


    ## Frame
    if (!exists("zero_vec")) {
      zero_vec <- 0
    }
    family_frame <- cbind(Y, mean_vec, theta_vec, zero_vec)
    if (y == "binomial") {
      pvec <- apply(family_frame, 1, function(x) {
        stats::pbinom(x[1], prob = x[2], size = 1)
      })
    } else if (y == "poisson") {
      pvec <- apply(family_frame, 1, function(x) {
        stats::ppois(x[1], lambda = x[2])
      })
    } else if (y == "gaussian") {
      pvec <- apply(family_frame, 1, function(x) {
        gamlss.dist::pNO(x[1], mu = x[2], sigma = abs(x[3]))
      })
    } else if (y == "nb") {
      pvec <- apply(family_frame, 1, function(x) {
        stats::pnbinom(x[1], mu = x[2], size = x[3])
      })
    } else if (y == "zip") {
      pvec <- apply(family_frame, 1, function(x) {
        gamlss.dist::pZIP(x[1], mu = x[2], sigma = abs(x[4]))
      })
    }
    else if (y == "zinb") {
      pvec <- apply(family_frame, 1, function(x) {
        gamlss.dist::pZINBI(x[1],
                            mu = x[2],
                            sigma = abs(x[3]),
                            nu = x[4])
      })
    } else {
      stop("Distribution of gamlss must be one of gaussian, binomial, poisson, nb, zip or zinb!")
    }

    ## CHECK ABOUT THE FIRST PARAM!!!!!
    if (DT) {
      if (y == "poisson") {
        pvec2 <- apply(family_frame, 1, function(x) {
          stats::ppois(x[1] - 1, lambda = x[2]) * as.integer(x[1] > 0)
        })
      } else if (y == "gaussian" | y == "binomial") {
        ## Gaussian is continuous, thus do not need DT.
        ## Binomial seems to be weird to have DT.
        message("Continuous gaussian does not need DT.")
        pvec2 <- pvec
        # pvec2 <- apply(family_frame, 1, function(x){
        # pNO(x[1]-1, mu = x[2], sigma = abs(x[3]))*as.integer(x[1] > 0)
        #})
      } else if (y == "nb") {
        pvec2 <- apply(family_frame, 1, function(x) {
          stats::pnbinom(x[1] - 1, mu = x[2], size = x[3]) * as.integer(x[1] > 0)
        })
      } else if (y == "zip") {
        pvec2 <- apply(family_frame, 1, function(x) {
          ifelse(x[1] > 0, gamlss.dist::pZIP(x[1] - 1, mu = x[2], sigma = abs(x[4])), 0)
        })
      } else if (y == "zinb") {
        pvec2 <- apply(family_frame, 1, function(x) {
          ifelse(x[1] > 0,
                 gamlss.dist::pZINBI(
                   x[1] - 1,
                   mu = x[2],
                   sigma = abs(x[3]),
                   nu = x[4]
                 ),
                 0)
        })
      } else {
        stop("Distribution of gamlss must be one of gaussian, binomial, poisson, nb, zip or zinb!")
      }

      u1 <- pvec
      u2 <- pvec2

      v <- stats::runif(length(mean_vec))
      ## Random mapping
      r <- u1 * v + u2 * (1 - v)
    } else {
      r <- pvec
    }

    if(length(r) < dim(sce)[2]){
      new_r <- rep(1, dim(sce)[2])
      names(new_r) <- colnames(sce)
      new_r[names(r)] <- r
      r <- new_r
    }

    ## Avoid Inf
    idx_adjust <- which(1 - r < epsilon)
    r[idx_adjust] <- r[idx_adjust] - epsilon
    idx_adjust <- which(r < epsilon)
    r[idx_adjust] <- r[idx_adjust] + epsilon

    if (pseudo_obs) {
      r <- rvinecopulib::pseudo_obs(r)
    }

    stats::qnorm(
      r,
      mean = 0,
      sd = 1,
      lower.tail = TRUE,
      log.p = FALSE
    )
  }

  paraFunc <- parallel::mcmapply
  if(.Platform$OS.type == "windows"){
    BPPARAM <- BiocParallel::SnowParam()
    parallelization <- "bpmapply"
  }
  if(parallelization == "bpmapply"){
    paraFunc <- BiocParallel::bpmapply
  }
  if(parallelization == "pbmcmapply"){
    paraFunc <- pbmcapply::pbmcmapply
  }
  if(parallelization == "bpmapply"){
    if(class(BPPARAM)[1] != "SerialParam"){
      BPPARAM$workers <- n_cores
    }
    mat <- paraFunc(mat_function, x = seq_len(dim(sce)[1]), y = family_use, SIMPLIFY = TRUE, BPPARAM = BPPARAM)
  }else{
    mat <- paraFunc(mat_function, x = seq_len(dim(sce)[1]), y = family_use, SIMPLIFY = TRUE, mc.cores = n_cores)
  }
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
                      data,
                      DT = TRUE,
                      pseudo_obs = FALSE,
                      epsilon = 1e-6,
                      n_cores,
                      family_use,
                      parallelization,
                      BPPARAM) {
  ## Extract count matrix
  count_mat <-
      t(as.matrix(SummarizedExperiment::assay(sce, assay_use)))
  removed_cell_list <- lapply(marginal_list, function(x){x$removed_cell})
  marginal_list <- lapply(marginal_list, function(x){x$fit})


  # n cell
  ncell <- dim(count_mat)[1]

  mat_function <- function(x, y) {
    fit <- marginal_list[[x]]
    removed_cell <- removed_cell_list[[x]]
    if(length(removed_cell) > 0 && !any(is.na(removed_cell))){
      data<- data[-removed_cell,]
    }
    if (methods::is(fit, "gamlss")) {
      mean_vec <- stats::predict(fit, type = "response", what = "mu", data = data) #
      if (y == "poisson" | y == "binomial") {
        theta_vec <- rep(NA, length(mean_vec))
      } else if (y == "gaussian") {
        theta_vec <-
          stats::predict(fit, type = "response", what = "sigma", data = data) # called the theta_vec but actually used as sigma_vec for Gaussian
      } else if (y == "nb") {
        theta_vec <-
          1 / stats::predict(fit, type = "response", what = "sigma", data = data) #
        #theta_vec[theta_vec < 1e-3] <- 1e-3
      } else if (y == "zip") {
        theta_vec <- rep(NA, length(mean_vec))
        zero_vec <-
          stats::predict(fit, type = "response", what = "sigma", data = data)
      } else if (y == "zinb") {
        theta_vec <- stats::predict(fit, type = "response", what = "sigma", data = data)
        zero_vec <-
          stats::predict(fit, type = "response", what = "nu", data = data)
      } else {
        stop("Distribution of gamlss must be one of gaussian, binomial, poisson, nb, zip or zinb!")
      }
    } else {
      ## if input is from mgcv, update the family
      y <- stats::family(fit)$family[1]
      if (grepl("Negative Binomial", y)) {
        y <- "nb"
      }

      mean_vec <- stats::predict(fit, type = "response")
      if (y == "poisson" | y == "binomial") {
        theta_vec <- rep(NA, length(mean_vec))
      } else if (y == "gaussian") {
        theta_vec <- rep(sqrt(fit$sig2), length(mean_vec)) # called the theta_vec but actually used as sigma_vec for Gaussian
      } else if (y == "nb") {
        theta <- fit$family$getTheta(TRUE)
        theta_vec <- rep(theta, length(mean_vec))
      } else {
        stop("Distribution of mgcv must be one of gaussian, binomial, poisson or nb!")
      }
    }

    ## Mean for Each Cell


    Y <- count_mat[names(mean_vec), x]


    ## Frame
    if (!exists("zero_vec")) {
      zero_vec <- 0
    }
    family_frame <- cbind(Y, mean_vec, theta_vec, zero_vec)

    if (y == "binomial") {
      pvec <- apply(family_frame, 1, function(x) {
        stats::pbinom(x[1], prob = x[2], size = 1)
      })
    } else if (y == "poisson") {
      pvec <- apply(family_frame, 1, function(x) {
        stats::ppois(x[1], lambda = x[2])
      })
    } else if (y == "gaussian") {
      pvec <- apply(family_frame, 1, function(x) {
        gamlss.dist::pNO(x[1], mu = x[2], sigma = abs(x[3]))
      })
    } else if (y == "nb") {
      pvec <- apply(family_frame, 1, function(x) {
        stats::pnbinom(x[1], mu = x[2], size = x[3])
      })
    } else if (y == "zip") {
      pvec <- apply(family_frame, 1, function(x) {
        gamlss.dist::pZIP(x[1], mu = x[2], sigma = abs(x[4]))
      })
    }
    else if (y == "zinb") {
      pvec <- apply(family_frame, 1, function(x) {
        gamlss.dist::pZINBI(x[1],
                            mu = x[2],
                            sigma = abs(x[3]),
                            nu = x[4])
      })
    } else {
      stop("Distribution of gamlss must be one of gaussian, poisson, nb, zip or zinb!")
    }

    ## CHECK ABOUT THE FIRST PARAM!!!!!
    if (DT) {
      if (y == "poisson") {
        pvec2 <- apply(family_frame, 1, function(x) {
          stats::ppois(x[1] - 1, lambda = x[2]) * as.integer(x[1] > 0)
        })
      } else if (y == "gaussian" | y == "binomial") {
        ## Gaussian is continuous, thus do not need DT.
        message("Continuous gaussian doesnot need DT.")
        pvec2 <- pvec
        # pvec2 <- apply(family_frame, 1, function(x){
        # pNO(x[1]-1, mu = x[2], sigma = abs(x[3]))*as.integer(x[1] > 0)
        #})
      } else if (y == "nb") {
        pvec2 <- apply(family_frame, 1, function(x) {
          stats::pnbinom(x[1] - 1, mu = x[2], size = x[3]) * as.integer(x[1] > 0)
        })
      } else if (y == "zip") {
        pvec2 <- apply(family_frame, 1, function(x) {
          ifelse(x[1] > 0, gamlss.dist::pZIP(x[1] - 1, mu = x[2], sigma = abs(x[4])), 0)
        })
      } else if (y == "zinb") {
        pvec2 <- apply(family_frame, 1, function(x) {
          ifelse(x[1] > 0,
                 gamlss.dist::pZINBI(
                   x[1] - 1,
                   mu = x[2],
                   sigma = abs(x[3]),
                   nu = x[4]
                 ),
                 0)
        })
      } else {
        stop("Distribution of gamlss must be one of gaussian, binomial, poisson, nb, zip or zinb!")
      }

      u1 <- pvec
      u2 <- pvec2

      v <- stats::runif(length(mean_vec))
      ## Random mapping
      r <- u1 * v + u2 * (1 - v)
    } else {
      r <- pvec
    }

    if(length(r) < dim(sce)[2]){
      new_r <- rep(1, dim(sce)[2])
      names(new_r) <- colnames(sce)
      new_r[names(r)] <- r
      r <- new_r
    }

    ## Avoid Inf
    idx_adjust <- which(1 - r < epsilon)
    r[idx_adjust] <- r[idx_adjust] - epsilon
    idx_adjust <- which(r < epsilon)
    r[idx_adjust] <- r[idx_adjust] + epsilon

    if (pseudo_obs) {
      r <- rvinecopulib::pseudo_obs(r)
    }

    r
  }

  paraFunc <- parallel::mcmapply
  if(.Platform$OS.type == "windows"){
    BPPARAM <- BiocParallel::SnowParam()
    parallelization <- "bpmapply"
  }
  if(parallelization == "bpmapply"){
    paraFunc <- BiocParallel::bpmapply
  }
  if(parallelization == "pbmcmapply"){
    paraFunc <- pbmcapply::pbmcmapply
  }
  if(parallelization == "bpmapply"){
    if(class(BPPARAM)[1] != "SerialParam"){
      BPPARAM$workers <- n_cores
    }
    mat <- paraFunc(mat_function, x = seq_len(dim(sce)[1]), y = family_use, SIMPLIFY = TRUE, BPPARAM = BPPARAM)
  }else{
    mat <- paraFunc(mat_function, x = seq_len(dim(sce)[1]), y = family_use, SIMPLIFY = TRUE, mc.cores = n_cores)
  }
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

cal_bic <- function(norm.mat,
                    cor.mat,
                    ind) {
  n_obs <- dim(norm.mat)[1]
  if (ind) {
    copula.bic = 0
  } else {
    copula.nop <- (as.integer(sum(cor.mat != 0)) - dim(cor.mat)[1]) / 2

    copula.bic <-
      -2 * (sum(
        mvtnorm::dmvnorm(
          x = norm.mat,
          mean = rep(0, dim(cor.mat)[1]),
          sigma = cor.mat,
          log = TRUE
        )
      ) - sum(rowSums(stats::dnorm(norm.mat, log = TRUE)))) + log(n_obs) * copula.nop
  }

  copula.bic
}

## Similar to the cora function from "Rfast" but uses different functions to calculate column means and row sums.
corrlation <- function(x) {
  mat <- t(x) - matrixStats::colMeans2(x)
  mat <- mat / sqrt(matrixStats::rowSums2(mat^2))
  tcrossprod(mat)
}

### Sample MVN based on cor
sampleMVN <- function(n,
                      Sigma, n_cores = n_cores, fastmvn = fastmvn) {
  if(fastmvn) {
    mvnrv <- mvnfast::rmvn(n, mu = rep(0, dim(Sigma)[1]), sigma = Sigma, ncores = n_cores)
  } else {
    mvnrv <-
      rmvnorm(n, mean = rep(0, dim(Sigma)[1]), sigma = Sigma, checkSymmetry = FALSE, method="eigen")
  }
  mvnrvq <- apply(mvnrv, 2, stats::pnorm)

  return(mvnrvq)
}

## fix integer overflow issue when # of gene* # of cell is too larger in rnorm(n * ncol(sigma))
## This function comes from package Mvnorm.
rmvnorm <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                    method=c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, checkSymmetry = TRUE)
{
  if (checkSymmetry && !isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                                    check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != nrow(sigma))
    stop("mean and sigma have non-conforming size")

  method <- match.arg(method)

  R <- if(method == "eigen") {
    ev <- eigen(sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))){
      warning("sigma is numerically not positive semidefinite")
    }
    ## ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
    ## faster for large  nrow(sigma):
    t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  }
  else if(method == "svd"){
    s. <- svd(sigma)
    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))){
      warning("sigma is numerically not positive semidefinite")
    }
    t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
  }
  else if(method == "chol"){
    R <- chol(sigma, pivot = TRUE)
    R[, order(attr(R, "pivot"))]
  }

  retval <- matrix(stats::rnorm(as.double(n) * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*%  R
  retval <- sweep(retval, 2, mean, "+")
  colnames(retval) <- names(mean)
  retval
}
#' Fit the marginal models
#'
#' \code{fit_marginal} fits the per-feature regression models.
#'
#' The function takes the result from \code{\link{construct_data}} as the input,
#' and fit the regression models for each feature based on users' specification.
#'
#' @param data An object from \code{\link{construct_data}}.
#' @param predictor A string of the predictor for the gam/gamlss model. Default is gene. This is essentially just a name.
#' @param mu_formula A string of the mu parameter formula
#' @param sigma_formula A string of the sigma parameter formula
#' @param family_use A string or a vector of strings of the marginal distribution.
#' Must be one of 'binomial', 'poisson', 'nb', 'zip', 'zinb' or 'gaussian', which represent 'poisson distribution',
#' 'negative binomial distribution', 'zero-inflated poisson distribution', 'zero-inflated negative binomail distribution',
#' and 'gaussian distribution' respectively.
#' @param n_cores An integer. The number of cores to use.
#' @param usebam A logic variable. If use \code{\link[mgcv]{bam}} for acceleration.
#' @param parallelization A string indicating the specific parallelization function to use.
#' Must be one of 'mcmapply', 'bpmapply', or 'pbmcmapply', which corresponds to the parallelization function in the package
#' \code{parallel},\code{BiocParallel}, and \code{pbmcapply} respectively. The default value is 'mcmapply'.
#' @param BPPARAM A \code{MulticoreParam} object or NULL. When the parameter parallelization = 'mcmapply' or 'pbmcmapply',
#' this parameter must be NULL. When the parameter parallelization = 'bpmapply',  this parameter must be one of the
#' \code{MulticoreParam} object offered by the package 'BiocParallel. The default value is NULL.
#' @param trace A logic variable. If TRUE, the warning/error log and runtime for gam/gamlss
#' will be returned, FALSE otherwise. Default is FALSE.
#' @return A list of fitted regression models. The length is equal to the total feature number.
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
#'
#' @export fit_marginal
#'
fit_marginal <- function(data,
                         predictor = "gene", ## Fix this later.
                         mu_formula,
                         sigma_formula,
                         family_use,
                         n_cores,
                         usebam,
                         parallelization = "mcmapply",
                         BPPARAM = NULL,
                         trace = FALSE) {

  count_mat <-  data$count_mat
  dat_cov <- data$dat
  feature_names <- colnames(count_mat)

  ## Check family_use
  if(length(family_use) == 1) {
    family_use <- rep(family_use, length(feature_names))
  }
  if(length(family_use) != length(feature_names)) {
    stop("The family_use must be either a single string or a vector with the same length as all features!")
  }


  fit_model_func <- function(gene,
                             family_gene,
                             dat_use,
                             #mgcv_formula,
                             mu_formula,
                             sigma_formula,
                             predictor,
                             count_mat
  ) {

    ## formula
    mgcv_formula <-
      stats::formula(paste0(predictor, "~", mu_formula))

    ## If use the mgcv s() smoother
    mu_mgcvform <- grepl("^s\\(", mu_formula) | grepl("^te\\(", mu_formula)

    ## If use bam to fit marginal distribution
    usebam <- usebam & mu_mgcvform ## If no smoothing terms, no need to to use bam.
    if(usebam){
      fitfunc = mgcv::bam
    }else{
      fitfunc = mgcv::gam
    }

    if (mu_mgcvform) {
      if(usebam){
        mu_formula <-
          stats::formula(paste0(predictor, "~", "ba(~", mu_formula, ", method = 'fREML', gc.level = 0, discrete = TRUE)"))
      }else{
        mu_formula <-
          stats::formula(paste0(predictor, "~", "ga(~", mu_formula, ", method = 'REML')"))
      }
    }
    else {
      mu_formula <- stats::formula(paste0(predictor, "~", mu_formula))
    }

    sigma_mgcvform <- grepl("^s\\(", sigma_formula) | grepl("^te\\(", sigma_formula)
    if (sigma_mgcvform) {
      if(usebam){
        sigma_formula <-
          stats::formula(paste0("~", "ba(~", sigma_formula, ", method = 'fREML', gc.level = 0, discrete = TRUE)"))
      }else{
        sigma_formula <-
          stats::formula(paste0("~", "ga(~", sigma_formula, ", method = 'REML')"))
      }

    } else {
      sigma_formula <- stats::formula(paste0("~", sigma_formula))
    }


    ## Add gene expr
    dat_use$gene <- count_mat[, gene]

    ## For error/warning logging
    add_log <- function(function_name, type, message) {
      new_l <- logs
      new_log <- list(function_name = function_name,
                      type = type,
                      message =  message)
      new_l[[length(new_l) + 1]]  <- new_log
      logs <<- new_l
    }

    logs <- list()
    ## Don't fit marginal if gene only have two or less non-zero expression
    if(length(which(dat_use$gene < 1e-5)) > length(dat_use$gene) - 2){
      add_log("fit_marginal","warning", paste0(gene, "is expressed in too few cells."))
      return(list(fit = NA, warning = logs, time = c(NA,NA)))
    }
    all_covariates <- all.vars(mgcv_formula)[-1]
    dat_cova <- dat_use[, all_covariates]
    check_factor <- all(sapply(dat_cova,is.factor))
    if (length(all_covariates) > 0 & check_factor){
      remove_idx_list <- lapply(all_covariates, function(x){
        curr_x <- tapply(dat_use$gene, dat_use[,x], sum)
        zero_group <- which(curr_x==0)
        if(length(zero_group) == 0){
          return(list(idx = NA, changeFormula = FALSE))
        }else{
          type <- names(curr_x)[zero_group]
          if(length(type) == length(unique(dat_use[,x])) - 1){
            return(list(idx = NA, changeFormula = TRUE))
          }
          return(list(idx = which(dat_use[,x] %in% type), changeFormula = FALSE))
        }

      })
      names(remove_idx_list) <- all_covariates
      remove_idx <- lapply(remove_idx_list, function(x)x$idx)
      remove_cell <- unlist(remove_idx)
      if(all(is.na(remove_cell))){
        remove_cell <- NA
      }else{
        remove_cell <- unique(stats::na.omit(remove_cell))
      }
      if(length(remove_cell) > 0 && !any(is.na(remove_cell))){
        dat_use <- dat_use[-remove_cell,]
      }

      changeFormula <-  sapply(remove_idx_list, function(x)x$changeFormula)
      if(length(which(changeFormula)) > 0){
        changeVars <- names(which(changeFormula))
        formulaUpdate <- paste0(changeVars, collapse = "-")
        mgcv_formula <- stats::update.formula(mgcv_formula, stats::as.formula(paste0("~.-",formulaUpdate)))
        mu_formula <- stats::update.formula(mu_formula, stats::as.formula(paste0("~.-",formulaUpdate)))
        sigmaVars <- which(changeVars %in% as.character(sigma_formula))
        if(length(sigmaVars) > 0){
          formulaUpdate <- paste0(changeVars[sigmaVars], collapse = "-")
        }
        sigma_formula = stats::update.formula(sigma_formula, stats::as.formula(paste0("~.-",formulaUpdate)))
      }

    }else{
      remove_cell <- NA
    }

    time_list <- c(NA,NA)

    if (family_gene == "binomial") {
      mgcv.fit <- withCallingHandlers(
        tryCatch({
          start.time <- Sys.time()
          res <-fitfunc(formula = mgcv_formula, data = dat_use, family = "binomial", discrete = usebam)
          end.time <- Sys.time()
          time <- as.numeric(end.time - start.time)
          time_list[1] <- time
          res
        }, error=function(e) {
          add_log("gam","error", toString(e))
          NULL
        }), warning=function(w) {
          add_log("gam","warning", toString(w))

        })

      if (sigma_formula != "~1") {
        gamlss.fit <- withCallingHandlers(
          tryCatch({
            start.time = Sys.time()
            res <- gamlss::gamlss(
              formula = mu_formula,
              #sigma.formula = sigma_formula, ## Binomial is one para dist.
              data = dat_use,
              family = gamlss.dist::BI,
              control = gamlss::gamlss.control(trace = FALSE, c.crit = 0.1)
            )
            end.time = Sys.time()
            time = as.numeric(end.time - start.time)
            time_list[2] <- time
            res
          }, error=function(e) {
            add_log("gamlss","error", toString(e))
            NULL
          }), warning=function(w) {
            add_log("gamlss","warning", toString(w))
          })

      } else {
        gamlss.fit <- NULL
      }
    } else if (family_gene == "poisson") {
      mgcv.fit <- withCallingHandlers(
        tryCatch({
          start.time <- Sys.time()
          res <-fitfunc(formula = mgcv_formula, data = dat_use, family = "poisson", discrete = usebam)
          end.time <- Sys.time()
          time <- as.numeric(end.time - start.time)
          time_list[1] <- time
          res
        }, error=function(e) {
          add_log("gam","error", toString(e))
          NULL
        }), warning=function(w) {
          add_log("gam","warning", toString(w))
        })

      if (sigma_formula != "~1") {
        gamlss.fit <- withCallingHandlers(
          tryCatch({
            start.time = Sys.time()
            res <- gamlss::gamlss(
              formula = mu_formula,
              #sigma.formula = sigma_formula, ## Poisson has constant mean
              data = dat_use,
              family = gamlss.dist::PO,
              control = gamlss::gamlss.control(trace = FALSE, c.crit = 0.1))
            end.time = Sys.time()
            time = as.numeric(end.time - start.time)
            time_list[2] <- time
            res
          }, error=function(e) {
            add_log("gamlss","error", toString(e))
            NULL
          }), warning=function(w) {
            add_log("gamlss","warning", toString(w))
          })
      } else {
        gamlss.fit <- NULL
      }
    } else if (family_gene == "gaussian") {
      mgcv.fit <- withCallingHandlers(
        tryCatch({
          start.time <- Sys.time()
          res <- fitfunc(formula = mgcv_formula, data = dat_use, family = "gaussian", discrete = usebam)
          end.time <- Sys.time()
          time <- as.numeric(end.time - start.time)
          time_list[1] <- time
          res
        }, error=function(e) {
          add_log("gam","error", toString(e))
          NULL
        }), warning=function(w) {
          add_log("gam","warning", toString(w))
        })

      if (sigma_formula != "~1") {
        gamlss.fit<- withCallingHandlers(
          tryCatch({
            start.time = Sys.time()
            res <- gamlss::gamlss(
              formula = mu_formula,
              sigma.formula = sigma_formula,
              data = dat_use,
              family = gamlss.dist::NO,
              control = gamlss::gamlss.control(trace = FALSE, c.crit = 0.1)
            )
            end.time = Sys.time()
            time = as.numeric(end.time - start.time)
            time_list[2] <- time
            res
          }, error=function(e) {
            add_log("gamlss","error", toString(e))
            NULL
          }), warning=function(w) {
            add_log("gamlss","warning", toString(w))
          })

      } else {
        gamlss.fit <- NULL
      }
    } else if (family_gene == "nb"){
      mgcv.fit <- withCallingHandlers(
        tryCatch({
          start.time <- Sys.time()
          res <-fitfunc(formula = mgcv_formula, data = dat_use, family = "nb", discrete = usebam)
          end.time <- Sys.time()
          time <- as.numeric(end.time - start.time)
          time_list[1] <- time
          res
        }, error=function(e) {
          add_log("gam","error", toString(e))
          NULL
        }), warning=function(w) {
          add_log("gam","warning", toString(w))
        })

      if (sigma_formula != "~1") {
        gamlss.fit <- withCallingHandlers(
          tryCatch({
            start.time = Sys.time()
            res <- gamlss::gamlss(
              formula = mu_formula,
              sigma.formula = sigma_formula,
              data = dat_use,
              family = gamlss.dist::NBI,
              control = gamlss::gamlss.control(trace = FALSE,  c.crit = 0.1)
            )
            end.time = Sys.time()
            time = as.numeric(end.time - start.time)
            time_list[2] <- time
            res
          }, error=function(e) {
            add_log("gamlss","error", toString(e))
            NULL
          }), warning=function(w) {
            add_log("gamlss","warning", toString(w))
          })
      } else {
        gamlss.fit <- NULL
      }
    } else if (family_gene == "zip") {
      mgcv.fit <- withCallingHandlers(
        tryCatch({
          start.time <- Sys.time()
          res <-fitfunc(formula = mgcv_formula, data = dat_use, family = "poisson", discrete = usebam)
          end.time <- Sys.time()
          time <- as.numeric(end.time - start.time)
          time_list[1] <- time
          res
        }, error=function(e) {
          add_log("gam","error", toString(e))
          NULL
        }), warning=function(w) {
          add_log("gam","warning", toString(w))
        })
      gamlss.fit <- withCallingHandlers(
        tryCatch({
          start.time = Sys.time()
          res <- gamlss::gamlss(
            formula = mu_formula,
            sigma.formula = mu_formula, ## Here sigma is the dropout prob, not variance!
            data = dat_use,
            family = gamlss.dist::ZIP,
            control = gamlss::gamlss.control(trace = FALSE, c.crit = 0.1)

          )
          end.time = Sys.time()
          time = as.numeric(end.time - start.time)
          time_list[2] <- time
          res
        }, error=function(e) {
          add_log("gamlss","error", toString(e))
          NULL
        }), warning=function(w) {
          add_log("gamlss","warning", toString(w))
        })

    } else if (family_gene == "zinb"){
      mgcv.fit <- withCallingHandlers(
        tryCatch({
          start.time <- Sys.time()
          res <- fitfunc(formula = mgcv_formula, data = dat_use, family = "nb", discrete = usebam)
          end.time <- Sys.time()
          time <- as.numeric(end.time - start.time)
          time_list[1] <- time
          res
        }, error=function(e) {
          add_log("gam","error", toString(e))
          NULL
        }), warning=function(w) {
          add_log("gam","warning", toString(w))
        })
      gamlss.fit <- withCallingHandlers(
        tryCatch({
          start.time = Sys.time()
          res <- gamlss::gamlss(
            formula = mu_formula,
            sigma.formula = sigma_formula,
            nu.formula = mu_formula, ## Here nu is the dropout probability!
            data = dat_use,
            family = gamlss.dist::ZINBI,
            control = gamlss::gamlss.control(trace = FALSE, c.crit = 0.1)
          )
          end.time = Sys.time()
          time = as.numeric(end.time - start.time)
          time_list[2] <- time
          res
        }, error=function(e) {
          add_log("gamlss","error", toString(e))
          NULL
        }), warning=function(w) {
          add_log("gamlss","warning", toString(w))
        })

    } else {
      stop("The regression distribution must be one of gaussian, poisson, nb, zip or zinb!")
    }

    ## Check if gamlss is fitted.
    if (!"gamlss" %in% class(gamlss.fit)) {
      if (sigma_formula != "~1") {
        message(paste0(gene, " uses mgcv::gam due to gamlss's error!"))
        ## gamlss.fit contains warning message
        if(!is.null(gamlss.fit)){
          ## check whether gam has warning messages
          if(is.null(warn)){
            warn = gamlss.fit
          }else{
            warn = c(warn, gamlss.fit)
          }
        }
      }

      fit <- mgcv.fit
    } else {

      mean_vec <- stats::predict(gamlss.fit, type = "response", what = "mu", data = dat_use)
      theta_vec <-
        stats::predict(gamlss.fit, type = "response", what = "sigma", data = dat_use)

      if_infinite <- (sum(is.infinite(mean_vec + theta_vec)) > 0)
      if_overmax <- (max(mean_vec, na.rm = TRUE) > 10* max(dat_use$gene, na.rm = TRUE))
      if(family_gene %in% c("nb","zinb")){
        #if_overdisp <- (min(theta_vec, na.rm = TRUE) < 1/ 1000)
        if_overdisp <- (max(theta_vec, na.rm = TRUE) > 1000)

      }else{
        if_overdisp <- FALSE
      }


      if (if_infinite | if_overmax | if_overdisp) {
        add_log("fit_marginal","warning", paste0(gene, " gamlss returns abnormal fitting values!"))
        #message(paste0(gene, " gamlss returns abnormal fitting values!"))
        fit <- mgcv.fit
      } else if (stats::AIC(mgcv.fit) - stats::AIC(gamlss.fit) < -Inf) {
        message(paste0(
          gene,
          "'s gamlss AIC is not signifincantly smaller than gam!"
        ))
        fit <- mgcv.fit
      }
      else {
        fit <- gamlss.fit
      }
    }

    if(trace){
      return(list(fit = fit, warning = logs, time = time_list, removed_cell = remove_cell))
    }
    return(list(fit = fit,removed_cell = remove_cell))
    #return(fit)
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
    model_fit <- suppressMessages(paraFunc(fit_model_func, gene = feature_names,
                                           family_gene = family_use,
                                           MoreArgs = list(dat_use = dat_cov,
                                                           #mgcv_formula = mgcv_formula,
                                                           mu_formula = mu_formula,
                                                           sigma_formula = sigma_formula,
                                                           predictor = predictor,
                                                           count_mat = count_mat),
                                           SIMPLIFY = FALSE, BPPARAM = BPPARAM))
  }else{
    model_fit <-  suppressMessages(paraFunc(fit_model_func, gene = feature_names,
                                           family_gene = family_use,
                                           mc.cores = n_cores,
                                           MoreArgs = list(dat_use = dat_cov,
                                                           #mgcv_formula = mgcv_formula,
                                                           mu_formula = mu_formula,
                                                           sigma_formula = sigma_formula,
                                                           predictor = predictor,
                                                           count_mat = count_mat),
                                           SIMPLIFY = FALSE))
  }

  # if(!is.null(model_fit$warning)) {
  #   #stop("Model has warning!")
  #   model_fit <- model_fit$value
  # }
  return(model_fit)
}

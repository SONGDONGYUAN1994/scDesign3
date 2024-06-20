#' Fit the marginal models
#'
#' \code{fit_marginal} fits the per-feature regression models.
#'
#' The function takes the result from \code{\link{construct_data}} as the input,
#' and fit the regression models for each feature based on users' specification.
#'
#' @param data An object from \code{\link{construct_data}}.
#' @param predictor A string of the predictor for the gam/gamlss model. Default is "gene". This is just a name.
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
#' @param trace A logic variable. If TRUE, the warning/error log and runtime for gam/gamlss will be returned.
#' will be returned, FALSE otherwise. Default is FALSE.
#' @param simplify A logic variable. If TRUE, the fitted regression model will only keep the essential contains for \code{predict}. Default is FALSE.
#' @param edf_flexible A logic variable. If TRUE, the fitted regression model will use the fitted relationship between gini coefficient and the effective degrees of freedom on a random selected gene sets. Default is FALSE.
#' @param seed A seed to control the random selection of genes.
#' 
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
                         predictor = "gene", 
                         mu_formula,
                         sigma_formula,
                         family_use,
                         n_cores,
                         usebam,
                         parallelization = "mcmapply",
                         BPPARAM = NULL,
                         trace = FALSE, 
                         simplify = FALSE,
                         edf_flexible = FALSE,  
                         seed = 123) {
  
  count_mat <-  data$count_mat
  dat_cov <- data$dat
  filtered_gene <- data$filtered_gene
  feature_names <- colnames(count_mat)
  
  
  # Extract K from mu formula
  matches <- regexpr("k\\s*=\\s*([0-9]+)", mu_formula, perl = TRUE)
  extracted_value <- regmatches(mu_formula, matches)
  extracted_K <- as.numeric(sub("k\\s*=\\s*", "", extracted_value))
  
  
  set.seed(seed)
  # Randomly select genes for edf fitting
  if(dim(count_mat)[2]>100 & extracted_K>200 & edf_flexible==TRUE){
    edf_fitting <- TRUE
    
    # genes for fitting edf-gini relationship
    edf_gini_genes <- sample(1:dim(count_mat)[2], 100)
    edf_gini_count_mat <-  count_mat[,edf_gini_genes]
    edf_gini_feature_names <- feature_names[edf_gini_genes]
    
    # genes for flexible edf
    edf_flexible_genes <- (1:dim(count_mat)[2])[-edf_gini_genes]
    edf_flexible_count_mat <- count_mat[,-edf_gini_genes]
    edf_flexible_feature_names <- feature_names[-edf_gini_genes]
    
  }else{
    edf_fitting <- FALSE
  }
  
  
  ## Check family_use
  if(length(family_use) == 1) {
    edf_gini_family_use <- rep(family_use, length(edf_gini_feature_names))
    edf_flexible_family_use <- rep(family_use, length(edf_flexible_feature_names))
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
                             count_mat,
                             edf=NULL
  ) {
    
    ## formula
    if(!is.null(edf)){
      mu_formula_ex <- sub("(k\\s*=).*", "\\1", mu_formula)
      mu_formula <- paste0("s(spatial1, spatial2, bs = 'gp', k= ",round(edf[[gene]]), ")")
    }
    
    mgcv_formula <-
      stats::formula(paste0(predictor, "~", mu_formula))
    
    ## If use the mgcv s() smoother
    mu_mgcvform <- grepl("s\\(", mu_formula) | grepl("te\\(", mu_formula)
    
    ## If use bam to fit marginal distribution
    usebam <- usebam & mu_mgcvform ## If no smoothing terms, no need to to use bam.
    if(usebam){
      fitfunc = mgcv::bam
    }else{
      fitfunc = mgcv::gam
    }
    
    if (mu_mgcvform) {
      terms <- attr(stats::terms(mgcv_formula), "term.labels")
      terms_smooth <- terms[which(grepl("s\\(", terms))] 
      
      if(usebam){
        terms_smooth_update <- sapply(terms_smooth, function(x){paste0("ba(~", x, ", method = 'fREML', gc.level = 0, discrete = TRUE)")})
        if(length(terms_smooth) == length(terms)){## only contain smooth terms
          mu_formula <-
            stats::formula(paste0(predictor, "~", paste0(terms_smooth_update, collapse = "+")))
        }else{
          terms_linear <- terms[which(!grepl("s\\(", terms))] 
          terms_update <- c(terms_linear, terms_smooth_update)
          mu_formula <-
            stats::formula(paste0(predictor, "~", paste0(terms_update, collapse = "+")))
        }
      }else{
        terms_smooth_update <- sapply(terms_smooth, function(x){paste0("ga(~", x, ", method = 'REML')")})
        if(length(terms_smooth) == length(terms)){## only contain smooth terms
          mu_formula <-
            stats::formula(paste0(predictor, "~", paste0(terms_smooth_update, collapse = "+")))
        }else{
          terms_linear <- terms[which(!grepl("s\\(", terms))] 
          terms_update <- c(terms_linear, terms_smooth_update)
          mu_formula <-
            stats::formula(paste0(predictor, "~", paste0(terms_update, collapse = "+")))
        }
      }
    }
    else {
      mu_formula <- stats::formula(paste0(predictor, "~", mu_formula))
    }
    
    sigma_mgcvform <- grepl("s\\(", sigma_formula) | grepl("te\\(", sigma_formula)
    if (sigma_mgcvform) {
      temp_sigma_formula <- stats::formula(paste0(predictor, "~", sigma_formula))
      terms <- attr(stats::terms(temp_sigma_formula), "term.labels")
      terms_smooth <- terms[which(grepl("s\\(", terms))]
      if(usebam){
        terms_smooth_update <- sapply(terms_smooth, function(x){paste0("ba(~", x, ", method = 'fREML', gc.level = 0, discrete = TRUE)")})
        if(length(terms_smooth) == length(terms)){## only contain smooth terms
          sigma_formula <-
            stats::formula(paste0("~", paste0(terms_smooth_update, collapse = "+")))
        }else{
          terms_linear <- terms[which(!grepl("s\\(", terms))] 
          terms_update <- c(terms_linear, terms_smooth_update)
          sigma_formula <-
            stats::formula(paste0("~", paste0(terms_update, collapse = "+")))
        }
      }else{
        terms_smooth_update <- sapply(terms_smooth, function(x){paste0("ga(~", x, ", method = 'REML')")})
        if(length(terms_smooth) == length(terms)){## only contain smooth terms
          sigma_formula <-
            stats::formula(paste0("~", paste0(terms_smooth_update, collapse = "+")))
        }else{
          terms_linear <- terms[which(!grepl("s\\(", terms))] 
          terms_update <- c(terms_linear, terms_smooth_update)
          sigma_formula <-
            stats::formula(paste0("~", paste0(terms_update, collapse = "+")))
        }
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
    if(!is.null(filtered_gene) & gene %in% filtered_gene){
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
    
    if(simplify) {
      fit <- simplify_fit(fit)
    }
    
    if(trace){
      return(list(fit = fit, warning = logs, time = time_list, removed_cell = remove_cell))
    }
    return(list(fit = fit,removed_cell = remove_cell))
    #return(fit)
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
  
  
  # If not using edf flexible fitting
  if(edf_fitting==FALSE){
    if(parallelization == "bpmapply"){
      if(class(BPPARAM)[1] != "SerialParam"){
        BPPARAM$workers <- n_cores
      }
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
  }else{ 
    # If using edf flexible fitting
    
    if(parallelization == "bpmapply"){
      if(class(BPPARAM)[1] != "SerialParam"){
        BPPARAM$workers <- n_cores
      }
      # Fit model to selected edf_gini_genes 
      model_fit_edf_gini <- suppressMessages(paraFunc(fit_model_func, gene = edf_gini_feature_names,
                                                      family_gene = edf_gini_family_use,
                                                      MoreArgs = list(dat_use = dat_cov,
                                                                      #mgcv_formula = mgcv_formula,
                                                                      mu_formula = mu_formula,
                                                                      sigma_formula = sigma_formula,
                                                                      predictor = predictor,
                                                                      count_mat = edf_gini_count_mat),
                                                      SIMPLIFY = FALSE, BPPARAM = BPPARAM))
    }else{
      
      # Fit model to selected edf_gini_genes
      model_fit_edf_gini <-  suppressMessages(paraFunc(fit_model_func, gene = edf_gini_feature_names,
                                                       family_gene = edf_gini_family_use,
                                                       mc.cores = n_cores,
                                                       MoreArgs = list(dat_use = dat_cov,
                                                                       #mgcv_formula = mgcv_formula,
                                                                       mu_formula = mu_formula,
                                                                       sigma_formula = sigma_formula,
                                                                       predictor = predictor,
                                                                       count_mat = edf_gini_count_mat),
                                                       SIMPLIFY = FALSE))
    }
    
    
    # Extract the fitted edf
    edf <- rep(NA, length(model_fit_edf_gini))
    for(i in 1:length(model_fit_edf_gini)){
      res_ind <- model_fit_edf_gini[i]
      if(lengths(res_ind)==2){
        res_ind <- res_ind[[names(res_ind)]]
        edf[i] <- sum(res_ind$fit$edf)
      }
    }
    
    # Fit a edf-gini relationship for edf_gini_genes
    edf_gini_count_gini <- apply(log(edf_gini_count_mat+1), MARGIN=2, FUN=reldist::gini)
    edf_gini_df <- data.frame(edf=edf, gini=edf_gini_count_gini)
    lm_edf_gini <- lm(edf~gini, data=edf_gini_df)
    # Upper bound for the lm coef
    #coef <- confint(lm_edf_gini)[,2]
    
    
    # Predict edf for edf_flexible_genes
    edf_flexible_count_gini <- apply(log(edf_flexible_count_mat+1), MARGIN=2, FUN=reldist::gini)
    edf_flexible_df <- data.frame(gini=edf_flexible_count_gini)
    edf_flexible_predicted <- predict(lm_edf_gini, edf_flexible_df, se.fit = TRUE, interval = "confidence", level = 0.95)
    edf_flexible_predicted_upr <- edf_flexible_predicted$fit[,3]
    
    
    # Fit again for the rest genes
    if(parallelization == "bpmapply"){
      if(class(BPPARAM)[1] != "SerialParam"){
        BPPARAM$workers <- n_cores
      }
      model_fit_edf_flexible <- suppressMessages(paraFunc(fit_model_func, gene = edf_flexible_feature_names,
                                                          family_gene = edf_flexible_family_use,
                                                          MoreArgs = list(dat_use = dat_cov,
                                                                          #mgcv_formula = mgcv_formula,
                                                                          mu_formula = mu_formula,
                                                                          sigma_formula = sigma_formula,
                                                                          predictor = predictor,
                                                                          count_mat = edf_flexible_count_mat,
                                                                          edf=edf_flexible_predicted_upr),
                                                          SIMPLIFY = FALSE, BPPARAM = BPPARAM))
    }else{
      model_fit_edf_flexible <-  suppressMessages(paraFunc(fit_model_func, gene = edf_flexible_feature_names,
                                                           family_gene = edf_flexible_family_use,
                                                           mc.cores = n_cores,
                                                           MoreArgs = list(dat_use = dat_cov,
                                                                           #mgcv_formula = mgcv_formula,
                                                                           mu_formula = mu_formula,
                                                                           sigma_formula = sigma_formula,
                                                                           predictor = predictor,
                                                                           count_mat = edf_flexible_count_mat,
                                                                           edf=edf_flexible_predicted_upr),
                                                           SIMPLIFY = FALSE))
    }
    
    
    # Combine model_fit_edf_gini and model_fit_edf_flexible
    model_fit <- vector(mode = "list", length = length(feature_names))
    names(model_fit) <- feature_names
    
    # Populate the new list based on indices:
    for (index in names(model_fit_edf_gini)) {
      model_fit[[index]] <- model_fit_edf_gini[[index]]
    }
    for (index in names(model_fit_edf_flexible)) {
      model_fit[[index]] <- model_fit_edf_flexible[[index]]
    }
  } 
  
  # if(!is.null(model_fit$warning)) {
  #   #stop("Model has warning!")
  #   model_fit <- model_fit$value
  # }
  return(model_fit)
}

simplify_fit <- function(cm) {
  ## This function is modified from https://win-vector.com/2014/05/30/trimming-the-fat-from-glm-models-in-r/
  cm$y = c()
  #cm$model = c()
  
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()  
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  
  #cm$mu.x = c()
  #cm$sigma.x = c()
  #cm$nu.x = c()
  
  #cm$family$variance = c()
  #cm$family$dev.resids = c()
  #cm$family$aic = c()
  #cm$family$validmu = c()
  #cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  
  attr(cm$mu.terms,".Environment") = c()
  attr(cm$mu.formula,".Environment") = c()
  
  attr(cm$sigma.terms,".Environment") = c()
  attr(cm$sigma.formula,".Environment") = c()
  
  attr(cm$nu.terms,".Environment") = c()
  attr(cm$nu.formula,".Environment") = c()
  cm
}

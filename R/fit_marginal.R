#' Fit the marginal models
#'
#' \code{fit_marginal} fits the per-feature regression models.
#'
#' The function takes the result from \code{\link{construct_data}} as the input,
#' and fit the regression models for each feature based on users' specification.
#'
#' @param data An object from \code{\link{construct_data}}. ## Fix later
#' @param predictor Default is gene. ## Fix later
#' @param mu_formula A string of the mu parameter formula
#' @param sigma_formula A string of the sigma parameter formula
#' @param family A string of the marginal distribution.
#' Must be one of 'poisson', 'nb', 'zip', 'zinb' or 'gaussian'.
#' @param n_cores An integer. The number of cores to use.
#' @param usebam A logic variable. If use \code{\link[mgcv]{bam}} for acceleration.
#'
#' @return A list of fitted regression models. The length is equal to the total feature number.
#'
#' @export fit_marginal
#'
fit_marginal <- function(data,
                         predictor = "gene", ## Fix this later.
                         mu_formula,
                         sigma_formula,
                         family,
                         n_cores,
                         usebam) {

  count_mat <-  data$count_mat
  dat <- data$dat
  feature_names <- colnames(count_mat)

  ## If use bam to fit marginal distribution
  if(usebam){
    fitfunc = mgcv::bam
  }else{
    fitfunc = mgcv::gam
  }


  mgcv_formula <-
    stats::as.formula(paste0(predictor, "~", mu_formula))


  ## If use the mgcv s() smoother
  mu_mgcvform <- grepl("^s\\(", mu_formula) | grepl("^te\\(", mu_formula)
  if (mu_mgcvform) {
    if(identical(fitfunc, mgcv::bam)){
      mu_formula <-
        stats::as.formula(paste0(predictor, "~", "ba(~", mu_formula, ", method = 'fREML', gc.level = 0, discrete = TRUE)"))
    }else{
      mu_formula <-
        stats::as.formula(paste0(predictor, "~", "ga(~", mu_formula, ", method = 'REML')"))
    }
  }
  else {
    mu_formula <- stats::as.formula(paste0(predictor, "~", mu_formula))
  }

  sigma_mgcvform <- grepl("^s\\(", sigma_formula) | grepl("^te\\(", sigma_formula)
  if (sigma_mgcvform) {
    if(identical(fitfunc, mgcv::bam)){
      sigma_formula <-
        stats::as.formula(paste0("~", "ba(~", sigma_formula, ", method = 'fREML', gc.level = 0, discrete = TRUE)"))
    }else{
      sigma_formula <-
        stats::as.formula(paste0("~", "ga(~", sigma_formula, ", method = 'REML')"))
    }

  } else {
    sigma_formula <- stats::as.formula(paste0("~", sigma_formula))
  }
  #
  # print(mu_formula)
  model_fit <- pbmcapply::pbmclapply(feature_names, function(gene,
                                                             mc.cores,
                                                             #pseudotime,
                                                             #celltype,
                                                             #spatial,
                                                             dat,
                                                             mgcv_formula,
                                                             mu_formula,
                                                             sigma_formula,
                                                             predictor,
                                                             count_mat,
                                                             family) {
    ## Add gene expr
    dat$gene <- count_mat[, gene]

    if (family == "poisson") {
      mgcv.fit <- fitfunc(formula = mgcv_formula, data = dat, family = "poisson")

      ## If sigma_formula == ~1, gamlss degenerates into mgcv::gam
      if (sigma_formula != "~1") {
        gamlss.fit <- tryCatch({
          res <- gamlss::gamlss(
            formula = mu_formula,
            #sigma.formula = sigma_formula, ## Poisson has constant mean
            data = dat,
            family = gamlss.dist::PO,
            control = gamlss::gamlss.control(trace = FALSE, c.crit = 0.1)

          )
          #message(paste0(gene, " gamlss fit successes!"))
          res
        },
        error = function(error) {
          message(paste0(gene, " gamlss fit fails!"))
          return(NULL)
        }#,
        #warning = function(warning) {
        #  message(paste0(gene, " gamlss fit ends with warning!"))
        #  return(NULL)
        #}
        , silent = FALSE)
      } else {
        gamlss.fit <- NULL
      }
    } else if (family == "gaussian") {
      # dat$gene = log1p(dat$gene)
      ## !!! Poisson doesnot have gamlss since its sigma equals mean!
      mgcv.fit <- fitfunc(formula = mgcv_formula, data = dat, family = "gaussian")

      if (sigma_formula != "~1") {
        gamlss.fit <- tryCatch({
          res <- gamlss::gamlss(
            formula = mu_formula,
            sigma.formula = sigma_formula,
            data = dat,
            family = gamlss.dist::NO,
            control = gamlss::gamlss.control(trace = FALSE, c.crit = 0.1)
          )
          #message(paste0(gene, " gamlss fit successes!"))
          res
        },
        error = function(error) {
          message(paste0(gene, " gamlss fit fails!"))
          return(NULL)
        }#,
        #warning = function(warning) {
        #  message(paste0(gene, " gamlss fit ends with warning!"))
        #  return(NULL)
        #}

        , silent = FALSE)
      } else {
        gamlss.fit <- NULL
      }
    } else if (family == "nb"){
      mgcv.fit <- fitfunc(formula = mgcv_formula, data = dat, family = "nb", discrete = identical(fitfunc, mgcv::bam))
      #mgcv.fit <- mgcv::gam(formula = mgcv_formula, data = dat, family = "nb")
      if (sigma_formula != "~1") {
        #dat$pseudotime <- pseudotime
        gamlss.fit <-         tryCatch({
          res <- gamlss::gamlss(
            formula = mu_formula,
            sigma.formula = sigma_formula,
            data = dat,
            family = gamlss.dist::NBI,
            control = gamlss::gamlss.control(trace = FALSE,  c.crit = 0.1)

          )
          #message(paste0(gene, " gamlss fit successes!"))
          res
        },
        error = function(error) {
          message(paste0(gene, " gamlss fit fails!"))
          return(NULL)
        }#,
        #warning = function(warning) {
        #  message(paste0(gene, " gamlss fit ends with warning!"))
        #  return(NULL)
        #}

        , silent = FALSE)
      } else {
        gamlss.fit <- NULL
      }
    } else if (family == "zip") {

      ## Fit mgcv::gam(poisson) in case gamlss(ZINB) fails.
      mgcv.fit <- fitfunc(formula = mgcv_formula, data = dat, family = "poisson")

      gamlss.fit <- tryCatch({
        res <- gamlss::gamlss(
          formula = mu_formula,
          sigma.formula = mu_formula, ## Here sigma is the dropout prob, not variance!
          data = dat,
          family = gamlss.dist::ZIP,
          control = gamlss::gamlss.control(trace = FALSE, c.crit = 0.1)

        )
        #message(paste0(gene, " gamlss fit successes!"))
        res
      },
      error = function(error) {
        message(paste0(gene, " gamlss fit fails!"))
        return(NULL)
      }#,
      #warning = function(warning) {
      #  message(paste0(gene, " gamlss fit ends with warning!"))
      #  return(NULL)
      #}

      , silent = FALSE)

    } else if (family == "zinb"){
      ## Fit mgcv::gam(poisson) in case gamlss(ZINB) fails.
      mgcv.fit <- fitfunc(formula = mgcv_formula, data = dat, family = "poisson")
      #mgcv.fit <- mgcv::gam(mgcv_formula, data = dat, family = "poisson")

      gamlss.fit <- tryCatch({
        res <- gamlss::gamlss(
          formula = mu_formula,
          sigma.formula = sigma_formula,
          nu.formula = mu_formula, ## Here nu is the dropout probability!
          data = dat,
          family = gamlss.dist::ZINBI,
          control = gamlss::gamlss.control(trace = FALSE, c.crit = 0.1)

        )
        #message(paste0(gene, " gamlss fit successes!"))
        res
      },
      error = function(error) {
        message(paste0(gene, " gamlss fit fails!"))
        return(NULL)
      }#,
      #warning = function(warning) {
      #  message(paste0(gene, " gamlss fit ends with warning!"))
      #  return(NULL)
      #}

      , silent = FALSE)
    } else {
      stop("The regression distribution must be one of gaussian, poisson, nb, zip or zinb!")
    }

    ## Check if gamlss is fitted.
    if (is.null(gamlss.fit)) {
      if ((sigma_formula != "~1")) {
        message(paste0(gene, " uses mgcv::gam due to gamlss's error!"))
      }

      fit <- mgcv.fit
    } else {
      mean_vec <- stats::predict(gamlss.fit, type = "response", what = "mu")
      theta_vec <-
        1 / stats::predict(gamlss.fit, type = "response", what = "sigma")

      if_infinite <- (sum(is.infinite(mean_vec + theta_vec)) > 0)
      if_overmax <- (max(mean_vec) > max(dat$gene))
      if (if_infinite | if_overmax) {
        message(paste0(gene, " gamlss returns abnormal fitting values!"))
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

    return(fit)
  },  mc.cores = n_cores,
  dat = dat,
  mgcv_formula = mgcv_formula,
  mu_formula = mu_formula,
  sigma_formula = sigma_formula,
  predictor = predictor,
  count_mat = count_mat,
  family = family)

  if(length(model_fit) == 2) {
    model_fit <- model_fit[[1]]
  }
  return(model_fit)
}

#' Perform the likelihood ratio test
#'
#' \code{perform_lrt} performs the likelihood ratio test to compare two list of marginal models.
#'
#' The function takes two lists of marginal models (by default, the first list is the alternative and the second is the null)
#' from \code{\link{fit_marginal}}. Note that LRT only makes sense for nested models.
#'
#' @param alter_marginal A list of marginal models from the alternative hypothesis.
#' @param null_marginal A list of marginal models from the null hypothesis. It must be strictly nested in the alternative model.
#'
#' @return A data.frame of the LRT result.
#'
#' @export perform_lrt

perform_lrt <- function(alter_marginal,
                        null_marginal) {

  p_tab <- mapply(function(fit1, fit2) {
    ## Check if same models for fit1/2
    sm <- identical(class(fit1), class(fit2))

    ## get predicted probabilities for both models
    m1y <- fit1$y
    m2y <- fit2$y
    m1n <- length(m1y)
    m2n <- length(m2y)
    if(m1n==0 | m2n==0)
      stop("Could not extract dependent variables from models.")

    if(m1n != m2n)
      stop(paste("Models appear to have different numbers of observations.\n",
                 "Model 1 has ",m1n," observations.\n",
                 "Model 2 has ",m2n," observations.\n",
                 sep="")
      )

    if(any(m1y != m2y)){
      stop(paste("Models appear to have different values on dependent variables.\n"))
    }

    ## gather up degrees of freedom
    if(methods::is(fit1, "gamlss")) {
      k1 <- fit1$df.fit
    } else if (methods::is(fit1, "gam"))
    {
      k1 <- sum(fit1$edf2 + fit1$edf1 - fit1$edf)
    } else {
      stop("Model must be either gamlss or mgcv::gam!")
    }

    if(methods::is(fit2, "gamlss")) {
      k2 <- fit2$df.fit
    } else if (methods::is(fit2, "gam"))
    {
      k2 <- sum(fit2$edf2 + fit2$edf1 - fit2$edf)
    } else {
      stop("Model must be either gamlss or mgcv::gam!")
    }

    #k1 <- length(coef(m1))
    #k2 <- length(coef(m2))

    ll1 <- stats::logLik(fit1)
    ll2 <- stats::logLik(fit2)

    sign_k <- sign(k1 - k2)
    lr <- -2*(ll2 - ll1)

    llk1 <- ll1/k1
    llk2 <- ll1/k2

    lrt.p <- ifelse(lr > 0, stats::pchisq(lr, (k1-k2), lower.tail = FALSE), NA)

    res <- c(sm, ll1, ll2, k1, k2, lrt.p)
    names(res) <- c("same_model", "LogLik_alter", "LogLik_null", "df_alter", "df_null", "p_value")
    return(res)
  }, fit1 = alter_marginal, fit2 = null_marginal)
  p_tab <- as.data.frame(t(p_tab))

  return(p_tab)
}

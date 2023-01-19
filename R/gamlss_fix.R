#' Functions from gamlss/gamlss.add with bugs fixed
#'
#' An additive function to be used while fitting GAMLSS models. The interface for \code{gam()} in the \pkg{mgcv} package.
#' @section ga
#' @param formula A formula of the model.
#' @param envir The environment.
#' @param control The control of the model fitting.
#' @param ... Other arguments.
#' @noRd
ga <- function(formula, envir, control = ga.control(...), ...)
{
  #------------------------------------------
  # function starts here
  #------------------------------------------
  scall <- deparse(sys.call(), width.cutoff = 500L)
  if (!methods::is(formula, "formula"))
    stop("formula argument in ga() needs a formula starting with ~")
  # get where "gamlss" is in system call, it can be in gamlss() or predict.gamlss()
  rexpr <- grepl("gamlss",sys.calls())
  #rexpr <- grepl("stats::model.frame.default", sys.calls()) ##



  for (i in length(rexpr):1) {
    position <- i # get the position
    if (rexpr[i] == TRUE)
      break
  }

  gamlss.env <- sys.frame(position) #gamlss or predict.gamlss
  #gamlss.env <- environment(formula)

  ## get the data
  ## this has been modified on the 12-12-14 to make sure that
  ##  if model.frame.gamlss() is used as for example in term.plot() the
  ## function does not fail (It need to implemented to all smoother using formula?)
  if (sys.call(position)[1] == "predict.gamlss()") {
    # if stats::predict is used
    Data <- get("data", envir = gamlss.env)
  } else if (sys.call(position)[1] == "gamlss()") {
    # if gamlss() is used
    if (is.null(get("gamlsscall", envir = gamlss.env)$data)) {
      # if no data argument but the formula can be interpreted
      Data <- stats::model.frame(formula)
    } else {
      # data argument in gamlss
      Data <- get("gamlsscall", envir = gamlss.env)$data
    }
  } else {
    Data <- get("data", envir = gamlss.env)
  }
  Data <- data.frame(base::eval(substitute(Data)))
  #-------------------------------------------------
  # new Daniil and Vlasis
  # Initialize gam
  formula <-
    stats::as.formula(paste0("Y.var", deparse(formula, width.cutoff = 500L)))

  Data$Y.var <- rep(0, nrow(Data))
  G <- mgcv::gam(
    formula,
    data = Data,
    offset = control$offset,
    method = control$method,
    optimizer = control$optimizer,
    control = control$control,
    scale =  control$scale,
    select = control$select,
    knots = control$knots,
    sp = control$sp,
    min.sp = control$min.sp,
    H = control$H,
    gamma = control$gamma,
    paraPen = control$paraPen,
    in.out = control$in.out,
    drop.unused.levels = control$drop.unused.levels,
    drop.intercept = control$drop.intercept,
    discrete = control$discrete,
    G = NULL,
    fit = FALSE
  )
  #--------------------------------------------------
  xvar <- rep(0,  dim(Data)[1])
  attr(xvar, "formula")     <- formula
  attr(xvar, "control")     <- control
  attr(xvar, "gamlss.env") <- gamlss.env
  attr(xvar, "data")       <- as.data.frame(Data)
  attr(xvar, "call")       <-
    substitute(gamlss.ga(data[[scall]], z, w, ...))
  attr(xvar, "class")      <- "smooth"
  attr(xvar, "G")          <- G
  xvar
}

#' Functions from gamlss/gamlss.add with bugs fixed
#'
#' An additive function to be used while fitting GAMLSS models. The interface for \code{bam()} in the \pkg{mgcv} package.
#' @section ba
#' @param formula A formula of the model.
#' @param control The control of the model fitting.
#' @param ... Other arguments.
#' @noRd
ba <-function(formula, control = ba.control(...), ...)
{
  #------------------------------------------
  # function starts here
  #------------------------------------------
  scall <- Reduce(paste, deparse(sys.call(), width.cutoff = 500L))
  if (!methods::is(formula, "formula"))
    stop("formula argument in ba() needs a formula starting with ~")
  # get where "gamlss" is in system call, it can be in gamlss() or predict.gamlss()
  rexpr <- grepl("gamlss",sys.calls()) ##
  #rexpr <- grepl("fitModel", sys.calls())
  #rexpr <- grepl("stats::model.frame.default", sys.calls())

  for (i in length(rexpr):1) {
    position <- i # get the position
    if (rexpr[i]==TRUE) break
  }
  gamlss.env <- sys.frame(position) #gamlss or predict.gamlss
  ## get the data
  ## this has been modified on the 12-12-14 to make sure that
  ##  if model.frame.gamlss() is used as for example in term.plot() the
  ## function does not fail (It need to implemented to all smoother using formulea?)
  if (sys.call(position)[1]=="predict.gamlss()") { # if stats::predict is used
    Data <- get("data", envir=gamlss.env)
  } else if (sys.call(position)[1]=="gamlss()") { # if gamlss() is used
    if (is.null(get("gamlsscall", envir=gamlss.env)$data)) { # if no data argument but the formula can be interpreted
      Data <- stats::model.frame(formula)
    } else {# data argument in gamlss
      Data <- get("gamlsscall", envir=gamlss.env)$data
    }
  } else {
    Data <- get("data", envir=gamlss.env)
  }
  Data <- data.frame(eval(substitute(Data)))
  #-------------------------------------------------
  # new Daniil and Vlasis
  # Initialize bam
  formula <- stats::as.formula(paste0("Y.var", Reduce(paste,deparse(formula, width.cutoff = 500L)) ))
  Data$Y.var = rep(0, nrow(Data))
  #browser()
  G = mgcv::bam(formula,
                data = Data,
                offset = control$offset,
                method = control$method,
                control = control$control,
                select = control$select,
                scale = control$scale,
                gamma = control$gamma,
                knots = control$knots,
                sp = control$sp,
                min.sp = control$min.sp,
                paraPen = control$paraPen,
                chunk.size = control$chunk.size,
                rho = control$rho,
                AR.start = control$AR.start,
                discrete = control$discrete,
                cluster = control$cluster,
                nthreads = control$nthreads,
                gc.level = control$gc.level,
                use.chol = control$use.chol,
                samfrac = control$samfrac,
                coef = control$coef,
                drop.unused.levels = control$drop.unused.levels,
                drop.intercept = control$drop.intercept,
                G = NULL,
                fit = FALSE
  )
  ##bam(formula, family=gaussian(),
  ##      data=list()#, weights=NULL, subset=NULL,
  #    na.action=na.omit, offset=NULL#, method="fREML"#,control=list()#,
  #    select=FALSE#, scale=0#,gamma=1,knots=NULL,sp=NULL,min.sp=NULL,
  #    paraPen=NULL,chunk.size=10000,rho=0,AR.start=NULL,discrete=FALSE,
  #    cluster=NULL,nthreads=1,gc.level=1,use.chol=FALSE,samfrac=1,
  #    coef=NULL,drop.unused.levels=TRUE,G=NULL,fit=TRUE,drop.intercept=NULL,...)
  #
  #--------------------------------------------------
  xvar <- rep(0,  dim(Data)[1])
  attr(xvar,"formula")     <- formula
  attr(xvar,"control")     <- control
  attr(xvar, "gamlss.env") <- gamlss.env
  attr(xvar, "data")       <- as.data.frame(Data)
  attr(xvar, "call")       <- substitute(gamlss.ba(data[[scall]], z, w, ...))
  attr(xvar, "class")      <- "smooth"
  attr(xvar, "G")          <- G
  xvar
}

#' Functions from gamlss/gamlss.add with bugs fixed
#'
#' The control for \code{ba()}. From \code{gamlss.add::ba.control()} and \code{gamlss::bam()}.
#' @section ba.control
#'
#' @param offset The offset in the formula.
#' @param method The method argument in \code{bam()}.
#' @param control A list of fit control parameters to replace defaults returned by gam.control. Any control parameters not supplied stay at their default values.
#' @param select The \code{select} argument in \code{bam()}. Determine should selection penalties be added to the smooth effects, so that they can in principle be penalized out of the model.
#' @param scale For the scale parameter. If this is positive then it is taken as the known scale parameter. Negative signals that the scale parameter is unknown. 0 signals that the scale parameter is 1 for Poisson and binomial and unknown otherwise.
#' @param gamma The \code{gamma} argument in \code{bam()}. Increase above 1 to force smoother fits.
#' @param knots The \code{knots} argument in \code{bam()}. An optional list containing user specified knot values to be used for basis construction.
#' @param sp The \code{sp} argument in \code{bam()}. A vector of smoothing parameters can be provided here.
#' @param min.sp The \code{min.sp} argument in \code{bam()}. Lower bounds can be supplied for the smoothing parameters.
#' @param paraPen The \code{paraPen} argument in \code{bam()}. Optional list specifying any penalties to be applied to parametric model terms.
#' @param chunk.size The model matrix is created in chunks of this size, rather than ever being formed whole.
#' @param rho An AR1 error model can be used for the residuals (based on dataframe order), of Gaussian-identity link models. This is the AR1 correlation parameter.
#' @param AR.start Logical variable of same length as data, \code{TRUE} at first observation of an independent section of AR1 correlation.
#' @param discrete With \code{method="fREML"} it is possible to discretize covariates for storage and efficiency reasons. If \code{discrete} is \code{TRUE}, a number or a vector of numbers for each smoother term, then discretization happens. If numbers are supplied they give the number of discretization bins.
#' @param cluster \code{bam} can compute the computationally dominant QR decomposition in parallel using parLapply from the \code{parallel} package, if it is supplied with a cluster on which to do this (a cluster here can be some cores of a single machine).
#' @param nthreads Number of threads to use for non-cluster computation (e.g. combining results from cluster nodes).
#' @param gc.level To keep the memory footprint down, it can help to call the garbage collector often, but this takes a substatial amount of time. Setting this to zero means that garbage collection only happens when R decides it should. Setting to 2 gives frequent garbage collection. 1 is in between.
#' @param use.chol By default \code{bam} uses a very stable QR update approach to obtaining the QR decomposition of the model matrix. For well conditioned models an alternative accumulates the crossproduct of the model matrix and then finds its Choleski decomposition, at the end. This is somewhat more efficient, computationally.
#' @param samfrac For very large sample size Generalized additive models the number of iterations needed for the model fit can be reduced by first fitting a model to a random sample of the data, and using the results to supply starting values. This initial fit is run with sloppy convergence tolerances, so is typically very low cost. \code{samfrac} is the sampling fraction to use. 0.1 is often reasonable.
#' @param coef Initial values for model coefficients.
#' @param drop.unused.levels By default unused levels are dropped from factors before fitting. For some smooths involving factor variables you might want to turn this off.
#' @param drop.intercept Set to \code{TRUE} to force the model to really not have the a constant in the parametric model part, even with factor variables present.
#' @param ... Other arguments.
#'
#' @return A control object
#' @noRd
ba.control <- function(offset = NULL,
                       method = "fREML",
                       control = list(),
                       select = FALSE,
                       scale = 0,
                       gamma = 1,
                       knots = NULL,
                       sp = NULL,
                       min.sp = NULL,
                       paraPen = NULL,
                       chunk.size = 10000,
                       rho = 0,
                       AR.start = NULL,
                       discrete = TRUE,
                       cluster = NULL,
                       nthreads = 2,
                       gc.level = 1,
                       use.chol = FALSE,
                       samfrac = 1,
                       coef = NULL,
                       drop.unused.levels = TRUE,
                       drop.intercept = NULL,
                       ...)
{
  #gam()
  control <- mgcv::gam.control(...)
  #ga()
  list( offset=offset, method=method, control=control,  select=select,
        scale=scale, gamma=gamma, knots=knots, sp=sp, min.sp=min.sp,
        paraPen = paraPen, chunk.size = chunk.size, rho = rho,
        AR.start = AR.start,
        discrete=discrete,  cluster=cluster, nthreads=nthreads,
        gc.level=gc.level, use.chol=use.chol, samfrac=samfrac,
        coef= coef,
        drop.unused.levels = drop.unused.levels,
        drop.intercept=drop.intercept, ...)
}

#' Functions from gamlss/gamlss.add with bugs fixed
#'
#' The control for \code{ga()}. From \code{gamlss.add::ga.control()} and \code{gamlss::gam()}.
#' @section ga.control
#' @param offset The offset in the formula.
#' @param method The smoothing parameter estimation method.
#' @param optimizer An array specifying the numerical optimization method to use to optimize the smoothing parameter estimation criterion (given by \code{method})
#' @param control A list of fit control parameters to replace defaults returned by \code{gam.control}.
#' @param scale If this is positive then it is taken as the known scale parameter. Negative signals that the scale parameter is unknown. 0 signals that the scale parameter is 1 for Poisson and binomial and unknown otherwise.
#' @param select If this is \code{TRUE} then \code{gam()} can add an extra penalty to each term so that it can be penalized to zero.
#' @param knots This is an optional list containing user specified knot values to be used for basis construction.
#' @param sp A vector of smoothing parameters can be provided here.
#' @param min.sp Lower bounds can be supplied for the smoothing parameters.
#' @param H A user supplied fixed quadratic penalty on the parameters of the GAM can be supplied, with this as its coefficient matrix.
#' @param gamma Increase this beyond 1 to produce smoother models.
#' @param paraPen Optional list specifying any penalties to be applied to parametric model terms.
#' @param in.out Optional list for initializing outer iteration.
#' @param drop.unused.levels By default unused levels are dropped from factors before fitting. For some smooths involving factor variables you might want to turn this off.
#' @param drop.intercept Set to \code{TRUE} to force the model to really not have the a constant in the parametric model part, even with factor variables present. Can be vector when \code{formula} is a list.
#' @param discrete Experimental option for setting up models for use with discrete methods employed in \code{bam}.
#' @param ... Other arguments
#'
#' @return A control object
#' @noRd
ga.control <- function(offset = NULL,
                       method = "REML",
                       optimizer = c("outer","newton"),
                       control = list(),
                       scale = 0,
                       select = FALSE,
                       knots = NULL,
                       sp = NULL,
                       min.sp = NULL,
                       H = NULL,
                       gamma = 1,
                       paraPen = NULL,
                       in.out = NULL,
                       drop.unused.levels = TRUE,
                       drop.intercept = NULL,
                       discrete = FALSE,
                       ...)
{
  #gam()
  control <- mgcv::gam.control(...)
  #ga()
  list(offset=offset, method=method, optimizer=optimizer, control=control,
       scale= scale,
       select=select, knots=knots, sp=sp, min.sp=min.sp, H=H, gamma=gamma,
       paraPen=paraPen, in.out=in.out,   drop.unused.levels = drop.unused.levels,
       drop.intercept=drop.intercept, discrete = discrete, ...)
}





#' Functions from gamlss/gamlss.add with bugs fixed
#'
#' The gamlss versions of the generic function \code{model.frame}
#' @section model.frame.gamlss
#' @param formula A formula of the model.
#' @param what For which parameter to extract the model.frame, terms or model.frame.
#' @param parameter Equivalent to \code{what}.
#' @param ... Other arguments.
#'
#' @return a vector or matrix of predicted values.
#' @noRd
model.frame.gamlss <- function(formula, what = c("mu", "sigma", "nu", "tau"), parameter = NULL, ...)
{
  object <- formula
  dots <- list(...)
  what <- if (!is.null(parameter)) {
    match.arg(parameter, choices = c("mu", "sigma", "nu", "tau"))
  } else match.arg(what)
  Call <- object$call
  parform <- stats::formula(object, what)
  data <- if (!is.null(Call$data)) {
    # problem here, as Call$data is .
    #eval(Call$data)
    # instead, this would work:
    if(what == "mu") {
      eval(Call$data, environment(formula$mu.terms))
    }
    else if (what == "sigma") {
      eval(Call$data, environment(formula$sigma.terms))
    } else if (what == "nu") {
      eval(Call$data, environment(formula$nu.terms))
    } else if (what == "tau") {
      eval(Call$data, environment(formula$tau.terms))
    }
    # (there is no formula$terms, just mu.terms and sigma.terms)
  } else {
    environment(formula$terms)
  }
  Terms <- stats::terms(parform)
  mf <- stats::model.frame(
    Terms,
    data,
    xlev = object[[paste(what, "xlevels", sep = ".")]]
  )
  mf
}

##' Support for Function ga()
##'
##'This is support for the  smoother functions \code{ga()} intefaces for Simon Woood's \code{gam()} functions from package \pkg{mgcv}. It is not intended to be called directly by users. From \code{gamlss.add::gamlss.ga}.
##' @param x The explanatory variables
##' @param y Iterative y variable
##' @param w Iterative weights
##' @param xeval If xeval=TRUE then predicion is used
##' @param ... Other arguments
##' @noRd
gamlss.ga <-function(x, y, w, xeval = NULL, ...) {
  if (is.null(xeval))
  {#fitting
    #formula <- attr(x,"formula")
    #control <- as.list(attr(x, "control"))
    Y.var <- y
    W.var <- w
    G <- attr(x,"G")
    G$y <- Y.var
    G$w <- W.var
    G$mf$Y.var <- Y.var
    G$mf$`(weights)` <- W.var
    fit <-  mgcv::gam(G=G, fit=TRUE)
    df <- sum(ifelse(is.null(fit$edf2), yes = fit$edf, fit$edf2) + fit$edf1 - fit$edf)-1
    fv <- stats::fitted(fit)
    residuals <- y-fv
    list(fitted.values=fv, residuals=residuals,
         nl.df = df, lambda=fit$sp[1], #
         coefSmo = fit, var=NA)    # var=fv has to fixed
  } else { # predict
    gamlss.env <- as.environment(attr(x, "gamlss.env"))
    obj <- get("object", envir=gamlss.env ) # get the object from predict
    TT <- get("TT", envir=gamlss.env ) # get wich position is now
    SL <- get("smooth.labels", envir=gamlss.env) # all the labels of the smoother
    fit <- eval(parse(text=paste("obj$", get("what", envir=gamlss.env), ".coefSmo[[",as.character(match(TT,SL)), "]]", sep="")))
    OData <- attr(x,"data")
    ll <- dim(OData)[1]
    pred <- stats::predict(fit,newdata = OData[seq(length(y)+1,ll),])
  }
}


##' Support for Function ba()
##'
##'This is support for the  smoother functions \code{ba()} intefaces for Simon Woood's \code{bam()} functions from package \pkg{mgcv}. It is not intended to be called directly by users. From \code{gamlss.add::gamlss.ba}.
##' @param x The explanatory variables
##' @param y Iterative y variable
##' @param w Iterative weights
##' @param xeval If xeval=TRUE then predicion is used
##' @param ... Other arguments
##' @noRd
gamlss.ba <-function(x, y, w, xeval = NULL, ...) {
  if (is.null(xeval))
  {#fitting
    Y.var <- y
    W.var <- w
    G <- attr(x,"G")
    control = attr(x,"control")
    G$y <- Y.var
    G$w <- W.var
    G$mf$Y.var <- Y.var
    G$mf$`(weights)` <- W.var
    fit <-  mgcv::bam(G=G, fit=TRUE,
                      offset=control$offset, method=control$method,
                      control=control$control, select=control$select,
                      scale=control$scale, gamma=control$gamma,
                      knots=control$knots, sp=control$sp, min.sp=control$min.sp,
                      paraPen=control$paraPen, chunk.size=control$chunk.size,
                      rho=control$rho, AR.start=control$AR.start,
                      discrete=control$discrete,
                      cluster=control$cluster, nthreads=control$nthreads,
                      gc.level=control$gc.level, use.chol=control$use.chol,
                      samfrac=control$samfrac,
                      drop.unused.levels=control$bam$drop.unused.levels)
    df <- sum(ifelse(is.null(fit$edf2), yes = fit$edf, fit$edf2) + fit$edf1 - fit$edf)-1
    fv <- stats::fitted(fit)
    residuals <- y-fv
    list( fitted.values=fv, residuals=residuals,
          nl.df = df, lambda=fit$sp[1], #
          coefSmo = fit, var=NA)    # var=fv has to fixed
  } else { # predict
    gamlss.env <- as.environment(attr(x, "gamlss.env"))
    obj <- get("object", envir=gamlss.env ) # get the object from predict
    TT <- get("TT", envir=gamlss.env ) # get wich position is now
    SL <- get("smooth.labels", envir=gamlss.env) # all the labels of the smoother
    fit <- eval(parse(text=paste("obj$", get("what", envir=gamlss.env), ".coefSmo[[",as.character(match(TT,SL)), "]]", sep="")))
    OData <- attr(x,"data")
    ll <- dim(OData)[1]
    pred <- stats::predict(fit,newdata = OData[seq(length(y)+1,ll),])
  }
}

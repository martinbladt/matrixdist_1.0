#' Regression method for ph Class
#'
#' @param x An object of class \linkS4class{ph}.
#' @param y Vector or data.
#' @param X Model matrix (no intercept needed).
#' @param B0 Initial regression coefficients (optional).
#' @param weight Vector of weights.
#' @param rcen Vector of right-censored observations.
#' @param rcenweight Vector of weights for right-censored observations.
#' @param stepsEM Number of EM steps to be performed.
#' @param methods Methods to use for matrix exponential calculation: `RM`, `UNI`, or `PADE`.
#' @param rkstep Runge-Kutta step size (optional).
#' @param uni_epsilon Epsilon parameter for uniformization method.
#' @param optim_method Method to use in gradient optimization.
#' @param maxit Maximum number of iterations when optimizing g function.
#' @param reltol Relative tolerance when optimizing g function.
#' @param every Number of iterations between likelihood display updates.
#'
#' @return An object of class \linkS4class{sph}.
#'
#' @importFrom methods is new
#' @importFrom stats optim
#' @importFrom utils tail
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' obj <- iph(ph(structure = "general", dimension = 2), gfun = "weibull", gfun_pars = 2)
#' data <- sim(obj, n = 100)
#' X <- runif(100)
#' reg(x = obj, y = data, X = X, stepsEM = 10)
setMethod(
  "reg", c(x = "ph", y = "ANY"),
  function(x,
           y,
           weight = numeric(0),
           rcen = numeric(0),
           rcenweight = numeric(0),
           X = numeric(0),
           B0 = numeric(0),
           stepsEM = 1000,
           methods = c("RK", "UNI"),
           rkstep = NA,
           uni_epsilon = NA,
           optim_method = "BFGS",
           maxit = 50,
           reltol = 1e-8,
           every = 10) {
    control <- if (optim_method == "BFGS") {
      list(
        maxit = maxit,
        factr = reltol,
        fnscale = -1
      )
    } else {
      list(
        maxit = maxit,
        reltol = reltol,
        fnscale = -1
      )
    }
    X <- as.matrix(X)
    if (methods[2] == "RK") {
      stop("For second method, select UNI or PADE (ordering avoided)")
    }
    if (any(dim(X) == 0)) {
      stop("input covariate matrix X, or use fit method instead")
    }
    is_iph <- is(x, "iph")
    rightCensored <- is.vector(rcen)
    
    # Distinguish functions for Right censoring and Interval censoring
    if(rightCensored){
      EMstep <- eval(parse(text = paste("EMstep_", methods[1], sep = "")))
    }else if(is.matrix(rcen)){
      EMstep <- eval(parse(text = "EMstep_UNI_intervalCensoring"))
    }
    
    if (!all(c(y, rcen) > 0)) {
      stop("data should be positive")
    }
    if (!all(c(weight, rcenweight) >= 0)) {
      stop("weights should be non-negative")
    }
    if (!is_iph) {
      par_g <- numeric(0)
      inv_g <- function(par, x) x
      if(rightCensored){
        LL_base <- eval(parse(text = paste("logLikelihoodPH_", methods[2], "s", sep = "")))
      }else{
        LL_base <- eval(parse(text = "logLikelihoodPH_UNIs_intervalCensoring"))
      }
      LL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
        ex <- exp(X %*% theta)
        scale1 <- ex[1:length(obs)]
        scale2 <-  if(rightCensored){tail(ex, length(rcens))}else{tail(ex, nrow(rcens))}
        return(LL_base(h, alpha, S, obs, weight, rcens, rcweight, scale1, scale2))
      }
    } else if (is_iph) {
      name <- x@gfun$name
      if (name %in% c("loglogistic", "gev")) {
        stop("not yet available for multi-parameter transforms")
      }
      par_g <- x@gfun$pars
      inv_g <- x@gfun$inverse
      if(rightCensored){
        LL_base <- eval(parse(text = paste("logLikelihoodM", x@gfun$name, "_", methods[2], "s", sep = "")))
      }else{
        LL_base <- eval(parse(text = paste("logLikelihoodM", x@gfun$name, "_UNIs_intervalCensoring", sep = "")))
      }
      LL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
        beta <- theta[1]
        B <- theta[2:length(theta)]
        if (beta < 0) {
          return(NA)
        }
        ex <- exp(X %*% B)
        scale1 <- ex[1:length(obs)]
        scale2 <-  if(rightCensored){tail(ex, length(rcens))}else{tail(ex, nrow(rcens))}
        return(LL_base(h, alpha, S, beta, obs, weight, rcens, rcweight, scale1, scale2))
      }
    }
    
    p <- dim(X)[2]
    n1 <- length(y)
    n2 <- if(rightCensored){ length(rcen)}else{nrow(rcen)}
    ng <- length(par_g)
    
    if (length(weight) == 0) weight <- rep(1, n1)
    if (length(rcenweight) == 0) rcenweight <- rep(1, n2)
    
    ph_par <- x@pars
    alpha_fit <- clone_vector(ph_par$alpha)
    S_fit <- clone_matrix(ph_par$S)
    if (length(B0) == 0) {
      B_fit <- rep(0, p)
    } else {
      B_fit <- B0
    }
    
    rcenk <- rcen
    rcenweightk <- rcenweight
    track <- numeric(stepsEM)
    
    options(digits.secs = 4)
    cat(format(Sys.time(), format = "%H:%M:%OS"), ": EM started", sep = "")
    cat("\n", sep = "")
    
    for (k in 1:stepsEM) {
      prop <- exp(X %*% B_fit)
      
      trans <- prop[1:n1] * inv_g(par_g, y)
      trans_cens <- prop[(n1 + 1):(n1 + n2)] * inv_g(par_g, rcen)
      
      A <- data_aggregation(trans, weight)
      if (length(rcen)>0) {
        B <- data_aggregation(rcen, rcenweight)
        rcenk <- B$un_obs
        rcenweightk <- B$weights
      }
      
      epsilon1 <- switch(which(methods[1] == c("RK", "UNI", "PADE")),
                         if (!is.na(rkstep)) rkstep else default_step_length(S_fit),
                         if (!is.na(uni_epsilon)) uni_epsilon else 1e-4,
                         0
      )
      epsilon2 <- switch(which(methods[2] == c("RK", "UNI", "PADE")),
                         if (!is.na(rkstep)) rkstep else default_step_length(S_fit),
                         if (!is.na(uni_epsilon)) uni_epsilon else 1e-4,
                         0
      )
      EMstep(epsilon1, alpha_fit, S_fit, A$un_obs, A$weights, rcenk, rcenweightk)
      theta <- c(par_g, B_fit)
      opt <- suppressWarnings(optim(
        par = theta,
        fn = LL,
        h = epsilon2,
        alpha = alpha_fit,
        S = S_fit,
        obs = y,
        weight = weight,
        rcens = rcen,
        rcweight = rcenweight,
        X = X,
        hessian = (k == stepsEM),
        method = optim_method,
        control = control
      ))
      
      track[k] <- opt$value
      par_g <- head(opt$par, ng)
      B_fit <- tail(opt$par, p)
      
      if (k %% every == 0) {
        cat("\r", "iteration:", k,
            ", logLik:", opt$value,
            sep = " "
        )
      }
    }
    cat("\n", sep = "")
    x@pars$alpha <- alpha_fit
    x@pars$S <- S_fit
    x@fit <- list(
      logLik = opt$value,
      nobs = sum(A$weights),
      hessian = opt$hessian,
      logLikHist = track
    )
    s <- sph(x, type = "reg")
    s@gfun$pars <- par_g
    s@coefs$B <- B_fit
    
    cat("\n", format(Sys.time(), format = "%H:%M:%OS"), ": EM finalized", sep = "")
    cat("\n", sep = "")
    
    s
  }
)

#' Aggregate data to improve computation speed
#' 
#' @param y Observations. Either a vector or a matrix.
#' @param w Respective weights of observations
#' 
#' @return Returns a named list with unique observations and associated weights. If y is a vector then the unique observations are given, otherwise the unique rows are returned.
#' 
data_aggregation <- function(y, w) {
  if ((length(w) == 0) && (is.vector(y))) w <- rep(1, length(y))
  if ((length(w) == 0) && (is.matrix(y))) w <- rep(1, nrow(y))
  
  observations <- cbind(y, w)
  mat <- data.frame(observations)
  
  if(is.vector(y)){
    names(mat) <- c("obs", "weight")
    agg <- stats::aggregate(mat$weight,
                            by = list(un_obs = mat$obs),
                            FUN = sum)
  }else if(is.matrix(y)){
    names(mat) <- c("lower", "upper", "weight")
    agg <- stats::aggregate(mat$weight,
                            by = list(lower = mat$lower, upper = mat$upper),
                            FUN = sum)
    
    agg$un_obs <- as.matrix(agg[, c("lower", "upper")])
  }
  
  list(un_obs = agg$un_obs, weights = agg$x)
}


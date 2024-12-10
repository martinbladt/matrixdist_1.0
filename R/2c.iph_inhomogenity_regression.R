#' Regression method for ph Class
#'
#' @param x An object of class \linkS4class{ph}.
#' @param y Vector or data.
#' @param X Model matrix (no intercept needed).
#' @param B0 Initial regression coefficients (optional).
#' @param X2 Model matrix for the inhomogeneity parameter (no intercept needed).
#' @param prop_f Regression function for the intensity function.
#' @param inhom_f Regression function for the intensity function.
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
#' inhomogeneity_reg(x = obj, y = data, X = X, stepsEM = 10)
setMethod(
  "inhomogeneity_reg", c(x = "ph", y = "ANY"),
  function(x,
           y,
           weight = numeric(0),
           rcen = numeric(0),
           rcenweight = numeric(0),
           X = numeric(0),
           X2 = NULL,
           prop_f = NULL,
           inhom_f = NULL,
           B0 = numeric(0),
           stepsEM = 1000,
           methods = c("RK", "UNI"),
           rkstep = NA,
           uni_epsilon = NA,
           optim_method = "BFGS",
           maxit = 50,
           reltol = 1e-8,
           break_tol = 1e-8,
           break_n = stepsEM,
           every = 10){
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
    # use same covariates for both proportionality term and intensity parameter
    if(is.null(X2)){X2 <- X}
    X <- as.matrix(X)
    X2 <- as.matrix(X2)
    
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
    
    # Base proportionality and intensity functions
    if(is.null(prop_f)){
      prop_f <- function(theta = NULL, data = X){rep(1, nrow(data))}
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
      LL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X, X2) {
        B <- head(theta, ncol(X)) 
        ex <- prop_f(theta = B, data = X)
        
        scale1 <- ex[1:length(obs)]
        scale2 <-  if(rightCensored){tail(ex, length(rcens))}else{tail(ex, nrow(rcens))}
        return(LL_base(h, alpha, S, obs, weight, rcens, rcweight, scale1, scale2))
      }
    } else if (is_iph) {
      
      if(is.null(inhom_f)){
        inhom_f <- function(phx = x, theta = NULL, data = X2){rep(phx@gfun$pars, nrow(data))}
      }
      
      name <- x@gfun$name
      if (name %in% c("loglogistic", "gev")) {
        stop("not yet available for multi-parameter transforms")
      }
      par_g <- x@gfun$pars
      inv_g <- x@gfun$inverse
      
      LL_base <- eval(parse(text = "logLikelihood_UNIs_PI"))
      
      LL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X, X2, gfun_name) {
        
        B <- head(theta, ncol(X)) 
        gamma <- tail(theta, ncol(X2)) 
        
        ex <- prop_f(theta = B, data = X) 
        beta <- inhom_f(theta = gamma, data = X2)
        
        if (any(beta < 0)) {
          return(NA)
        }
        
        scale1 <- ex[1:length(obs)]
        scale2 <-  if(rightCensored){tail(ex, length(rcens))}else{tail(ex, nrow(rcens))}
        
        beta1 <- head(beta,length(obs)) 
        beta2 <- if(rightCensored){tail(beta, length(rcens))}else{tail(beta, nrow(rcens))}
        
        return(LL_base(h, alpha, S, beta1, beta2, obs, weight, rcens, rcweight, scale1, scale2, gfun_name))
      }
    }
    
    p0 <- dim(X)[2]
    p1 <- dim(X2)[2]
    p <- p0 + p1
    
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
      # proportionality and inhomogeneity terms
      prop <- prop_f(theta = head(B_fit, p0), data = X)
      beta_k <- inhom_f(theta = tail(B_fit, p1), data = X2)
      
      # transformed observations
      trans <- prop[1:n1] * inv_g(beta_k[1:n1], y)
      trans_cens <- prop[(n1 + 1):(n1 + n2)] * inv_g(beta_k[(n1 + 1):(n1 + n2)], rcen)
      
      A <- data_aggregation(trans, weight)
      if (length(rcen)>0) {
        Bcens <- data_aggregation(rcen, rcenweight)
        rcenk <- Bcens$un_obs
        rcenweightk <- Bcens$weights
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
      theta <- B_fit
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
        X2 = X2,
        gfun_name = name,
        hessian = (k == stepsEM),
        method = optim_method,
        control = control
      ))

      track[k] <- opt$value
      B_fit <- opt$par
      
      if (k %% every == 0) {
        cat("\r", "iteration:", k,
            ", logLik:", opt$value,
            sep = " "
        )
      }
      
      # stop the algorithm if it is stuck somewhere
      if(k > break_n){
        b_diff <- track[k]-track[k-1]
        if(b_diff < break_tol){break}
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
    
    s@gfun$pars <- inhom_f(theta = tail(B_fit, p1), data = X2)
    s@gfun$prop <- prop_f(theta = head(B_fit, p0), data = X)
    
    s@gfun$regression_intensity <- function(t, beta = s@gfun$pars, prop = s@gfun$prop){
      s@gfun$intensity(beta, t)*prop
    }
    
    s@coefs$B_prop <- head(B_fit, p0)
    s@coefs$B_inhom <- tail(B_fit, p1)
    
    cat("\n", format(Sys.time(), format = "%H:%M:%OS"), ": EM finalized", sep = "")
    cat("\n", sep = "")
    
    s
  })

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

#' Regression Method for ph Class
#'
#' @param x an object of class \linkS4class{ph}.
#' @param y vector or data.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#'
setMethod(
  "reg2", c(x = "iph", y = "ANY"),
  function(x,
           y,
           weight = numeric(0),
           rcen = numeric(0),
           rcenweight = numeric(0),
           X = numeric(0),
           Z = numeric(0),
           stepsEM = 1000,
           B0 = numeric(0),
           C0 = numeric(0)) {
    X <- as.matrix(X)
    Z <- as.matrix(X)
    if(any(dim(Z) == 0)) Z <- X
    if(any(dim(X) == 0)) stop("input covariate matrix X, or use fit method instead")
    name <- x@gfun$name
    par_g <- x@gfun$pars
    specs <- reg2_g_specs(name) 
    inv_g <- specs$inv_g 
    mLL <- specs$mLL
    
    p1 <- dim(X)[2]
    p2 <- dim(Z)[2]
    n1 <- length(y)
    n2 <- length(rcen)

    if(length(weight) == 0) weight <- rep(1, n1)
    if(length(rcenweight) == 0) rcenweight <- rep(1, n2)
    
    ph_par <- x@pars
    pi_fit <- clone_vector(ph_par$alpha)
    T_fit <- clone_matrix(ph_par$S)
    
    if(length(B0) == 0){B_fit <- rep(0, p1)
    }else{B_fit <- B0}
    if(length(C0) == 0){C_fit <- rep(0, p2)
    }else{C_fit <- C0}
    C_intercept <- log(par_g)
    
    for (k in 1:stepsEM) {
      prop <- exp(X%*%B_fit)
      par_g <- exp(C_intercept + Z%*%C_fit) 
      
      trans <- inv_g(y, weight, par_g[1:n1]); trans$obs <- prop[1:n1] * trans$obs
      trans_cens <- inv_g(rcen, rcenweight, par_g[(n1 + 1):(n1 + n2)]); trans_cens$obs <- prop[(n1 + 1):(n1 + n2)] * trans_cens$obs
      
      A <- data_aggregation(trans$obs, trans$weight)
      B <- data_aggregation(trans_cens$obs, trans_cens$weight)
      
      RKstep <- default_step_length(T_fit)
      EMstep_RK(RKstep, pi_fit, T_fit, A$un_obs, A$weights, B$un_obs, B$weights)
      theta <- c(B_fit, C_intercept, C_fit)
      opt <- suppressWarnings(optim(par = theta, fn = mLL, 
                                    h = RKstep, 
                                    alpha = pi_fit, 
                                    S = T_fit, 
                                    obs = y,
                                    weight = weight, 
                                    rcens = rcen, 
                                    rcweight = rcenweight,
                                    X = X,
                                    Z = Z,
                                    p1 = p1,
                                    p2 = p2,
                                    hessian = (k == stepsEM),
                                    method = ifelse(k == stepsEM, "Nelder-Mead", "Nelder-Mead"),
                                    control = list(
                                      maxit = ifelse(k == stepsEM, 1000, 50),
                                      reltol = ifelse(k == stepsEM, 1e-8, 1e-6)
                                    )
                                    ))
      B_fit <- head(opt$par, p1)
      C_intercept <- opt$par[p1 + 1]
      C_fit <- tail(opt$par, p2)
      if (k %% 10 == 0) {
        cat("\r", "iteration:", k,
            ", logLik:", - opt$value,
            sep = " ")
      }
    }
    cat("\n", sep = "")
    x@pars$alpha <- pi_fit
    x@pars$S <- T_fit
    x@gfun$pars <- exp(C_intercept)
    x@fit <- list(cov = safe_cov(opt$hessian),
                  loglik = - opt$value)
    s <- sph(x, type = "reg2")
    s@coefs$B <- B_fit
    s@coefs$C <- C_fit
    return(s)
  }
)

reg2_g_specs <- function(name){
   if(name == "Weibull"){
    inv_g <- function(t, w, beta) return(list(obs = t^{beta}, weight = w)) 
    mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X, Z, p1, p2) {
      B <- theta[1:p1]
      C <- tail(theta, p2 + 1)
      beta <- exp(C[1] + Z%*%C[-1])
      ex <- exp(X%*%B)
      scale1 <- ex[1:length(obs)]; beta1 <- beta[1:length(obs)]
      scale2 <- tail(ex, length(rcens)); beta2 <- tail(beta, length(rcens))
      o1 <- order(scale1 * inv_g(obs, weight, beta1)$obs)
      o2 <- order(scale2 * inv_g(rcens, rcweight, beta2)$obs)
      return(- logLikelihoodMWeib_RKs_double(h, alpha, S, beta1[o1], beta2[o2], obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
    }
  }
  else if(name == "Pareto"){
    inv_g <- function(t, w, beta) return(list(obs = log(t/beta + 1), weight = w))
    mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X, Z, p1, p2) {
      B <- theta[1:p1]
      C <- tail(theta, p2 + 1)
      beta <- exp(C[1] + Z%*%C[-1])
      ex <- exp(X%*%B)
      scale1 <- ex[1:length(obs)]; beta1 <- beta[1:length(obs)]
      scale2 <- tail(ex, length(rcens)); beta2 <- tail(beta, length(rcens))
      o1 <- order(scale1 * inv_g(obs, weight, beta1)$obs)
      o2 <- order(scale2 * inv_g(rcens, rcweight, beta2)$obs)
      return(- logLikelihoodMPar_RKs_double(h, alpha, S, beta1[o1], beta2[o2], obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
    }
  }
  else if(name == "Gompertz"){
    inv_g <- function(t, w, beta) return(list(obs = (exp(t * beta) - 1) / beta, weight = w))
    mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X, Z, p1, p2) {
      B <- theta[1:p1]
      C <- tail(theta, p2 + 1)
      beta <- exp(C[1] + Z%*%C[-1])
      ex <- exp(X%*%B)
      scale1 <- ex[1:length(obs)]; beta1 <- beta[1:length(obs)]
      scale2 <- tail(ex, length(rcens)); beta2 <- tail(beta, length(rcens))
      o1 <- order(scale1 * inv_g(obs, weight, beta1)$obs)
      o2 <- order(scale2 * inv_g(rcens, rcweight, beta2)$obs)
      return(- logLikelihoodMGomp_RKs_double(h, alpha, S, beta1[o1], beta2[o2], obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
    }
  }
  else if(name == "LogNormal"){
    inv_g <- function(t, w, beta) return(list(obs = log(t + 1)^{beta}, weight = w)) 
    mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X, Z, p1, p2) {
      B <- theta[1:p1]
      C <- tail(theta, p2 + 1)
      beta <- exp(C[1] + Z%*%C[-1])
      ex <- exp(X%*%B)
      scale1 <- ex[1:length(obs)]; beta1 <- beta[1:length(obs)]
      scale2 <- tail(ex, length(rcens)); beta2 <- tail(beta, length(rcens))
      o1 <- order(scale1 * inv_g(obs, weight, beta1)$obs)
      o2 <- order(scale2 * inv_g(rcens, rcweight, beta2)$obs)
      return(- logLikelihoodMLogNormal_RKs_double(h, alpha, S, beta1[o1], beta2[o2], obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
    }
  }
  else{
    stop("fit for this gfun is not yet implemented")
  }
  return(list(inv_g = inv_g, mLL = mLL))
}


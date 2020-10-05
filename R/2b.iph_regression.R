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
  "reg", c(x = "ph", y = "ANY"),
  function(x,
           y,
           weight = numeric(0),
           rcen = numeric(0),
           rcenweight = numeric(0),
           X = numeric(0),
           stepsEM = 1000) {
    X <- as.matrix(X)
    if(any(dim(X) == 0)) stop("input covariate matrix X, or use fit method instead")
    is_iph <- is(x, "iph")
    if(is_iph){
      name <- x@gfun$name
      par_g <- x@gfun$pars
      specs <- reg_g_specs(name) 
      inv_g <- specs$inv_g 
      mLL <- specs$mLL
    }else{
      par_g <- numeric(0)
      specs <- reg_g_specs("Homogeneous") 
      inv_g <- specs$inv_g 
      mLL <- specs$mLL
    }
    A <- data_aggregation(y, weight); y <- A$un_obs; weight <- A$weights
    B <- data_aggregation(rcen, rcenweight); rcen <- B$un_obs; rcenweight <- B$weights
    
    p <- dim(X)[2]
    n1 <- length(y)
    n2 <- length(rcen)
    ng <- length(par_g)

    ph_par <- x@pars
    pi_fit <- clone_vector(ph_par$alpha)
    T_fit <- clone_matrix(ph_par$S)
    B_fit <- rep(0, p)

    for (k in 1:stepsEM) {
      prop <- exp(X%*%B_fit)
      
      trans <- inv_g(y, weight, par_g); trans$obs <- prop[1:n1] * trans$obs
      o1 <- order(trans$obs)
      trans <- lapply(trans, function(x) x[o1])
      
      trans_cens <- inv_g(rcen, rcenweight, par_g); trans_cens$obs <- prop[(n1 + 1):(n1 + n2)] * trans_cens$obs
      o2 <- order(trans_cens$obs)
      trans_cens <- lapply(trans_cens, function(x) x[o2])
      
      RKstep <- default_step_length(T_fit)
      EMstep_RK(RKstep, pi_fit, T_fit, trans$obs, trans$weight, trans_cens$obs, trans_cens$weight)
      theta <- c(par_g, B_fit)
      opt <- suppressWarnings(optim(par = theta, fn = mLL, 
                                    h = RKstep, 
                                    alpha = pi_fit, 
                                    S = T_fit, 
                                    obs = y,
                                    weight = weight, 
                                    rcens = rcen, 
                                    rcweight = rcenweight,
                                    X = X))
      par_g <- head(opt$par, ng)
      B_fit <- tail(opt$par, p)
      if (k %% 10 == 0) {
        cat("\r", "iteration:", k,
            ", logLik:", - opt$value,
            sep = " ")
      }
    }
    cat("\n", sep = "")
    x@pars$alpha <- pi_fit
    x@pars$S <- T_fit
    if(is_iph) x <- iph(x, gfun = name, gfun_pars = par_g)
    cat("\n", "regression parameters:", "\n", sep = " ")
    print(B_fit)
    return(x)
  }
)

reg_g_specs <- function(name){
  if(name == "Homogeneous"){
    inv_g <- function(t, w, beta) return(list(obs = t, weight = w)) 
    mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
      ex <- exp(X%*%theta)
      scale1 <- ex[1:length(obs)]
      scale2 <- tail(ex, length(rcens))
      o1 <- order(scale1 * obs)
      o2 <- order(scale2 * rcens)
      return(- logLikelihoodPH_RKs(h, alpha, S, obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
    }
  }
  # if(name == "Homogeneous"){
  #   inv_g <- function(t, w, beta) return(list(obs = t, weight = w)) 
  #   mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
  #     ex <- exp(X%*%theta)
  #     scale1 <- ex[1:length(obs)]
  #     scale2 <- tail(ex, length(rcens))
  #     return(- logLikelihoodPH_RKs2(h, alpha, S, obs, weight, rcens, rcweight, scale1, scale2))
  #   }
  # }
  else if(name == "Weibull"){
    inv_g <- function(t, w, beta) return(list(obs = t^{beta}, weight = w)) 
    mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
      beta <- theta[1]; B <- theta[2:length(theta)]
      if(beta < 0) return(NA)
      ex <- exp(X%*%B)
      scale1 <- ex[1:length(obs)]
      scale2 <- tail(ex, length(rcens))
      o1 <- order(scale1 * obs)
      o2 <- order(scale2 * rcens)
      return(- logLikelihoodMWeib_RKs(h, alpha, S, beta, obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
    }
  }
  else if(name == "Pareto"){
    inv_g <- function(t, w, beta) return(list(obs = log(t/beta + 1), weight = w))
    mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
      beta <- theta[1]; B <- theta[2:length(theta)]
      if(beta < 0) return(NA)
      ex <- exp(X%*%B)
      scale1 <- ex[1:length(obs)]
      scale2 <- tail(ex, length(rcens))
      o1 <- order(scale1 * obs)
      o2 <- order(scale2 * rcens)
      return(- logLikelihoodMPar_RKs(h, alpha, S, beta, obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
    }
  }
  else if(name == "LogLogistic"){
    inv_g <- function(t, w, beta) return(list(obs = log((t/beta[1])^{beta[2]} + 1), weight = w))
    mLL <- function(h, alpha, S, beta, obs, weight, rcens, rcweight, X) {
      beta <- theta[1:2]; B <- theta[3:length(theta)]
      if(beta[1] < 0 | beta[2] < 0) return(NA)
      ex <- exp(X%*%B)
      scale1 <- ex[1:length(obs)]
      scale2 <- tail(ex, length(rcens))
      o1 <- order(scale1 * obs)
      o2 <- order(scale2 * rcens)
      return(- logLikelihoodMLogLogistic_RKs(h, alpha, S, beta, obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
    }
  }
  else if(name == "Gompertz"){
    inv_g <- function(t, w, beta) return(list(obs = (exp(t * beta) - 1) / beta, weight = w))
    mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
      beta <- theta[1]; B <- theta[2:length(theta)]
      if(beta < 0) return(NA)
      ex <- exp(X%*%B)
      scale1 <- ex[1:length(obs)]
      scale2 <- tail(ex, length(rcens))
      o1 <- order(scale1 * obs)
      o2 <- order(scale2 * rcens)
      return(- logLikelihoodMGomp_RKs(h, alpha, S, beta, obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
    }
  }
  else if(name == "GEVD"){
    inv_g <- reversTransformData
    mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
      beta <- theta[1:3]; B <- theta[4:length(theta)]
      if(beta[2] < 0) return(NA)
      ex <- exp(X%*%B)
      scale1 <- ex[1:length(obs)]
      scale2 <- tail(ex, length(rcens))
      o1 <- order(scale1 * obs)
      o2 <- order(scale2 * rcens)
      return(- logLikelihoodMGEV_RKs(h, alpha, S, beta, obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
    }
  }else{
    stop("fit for this gfun is not yet implemented")
  }
  return(list(inv_g = inv_g, mLL = mLL))
}


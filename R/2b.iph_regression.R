#' Regression Method for ph Class
#'
#' @param x an object of class \linkS4class{ph}.
#' @param y vector or data.
#' @param X model matrix (no intercept needed).
#' @param B0 initial regression coefficients (optional).
#' @param weight vector of weights.
#' @param rcen vector of right-censored observations
#' @param rcenweight vector of weights for right-censored observations.
#' @param stepsEM number of EM steps to be performed.
#' @param methods methods to use for matrix exponential calculation: RM, UNI or PADE
#' @param rkstep Runge-Kutta step size (optional)
#' @param uni_epsilon epsilon parameter for uniformization method
#' @param optim_method method to use in gradient optimization
#' @param maxit maximum number of iterations when optimizing g function.
#' @param reltol relative tolerance when optimizing g function.
#' @param every number of iterations between likelihood display updates.
#'
#' @return An object of class \linkS4class{sph}.
#' 
#' @importFrom methods is new
#' @importFrom stats optim
#' @importFrom utils tail
#' 
#' @export
#'
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
    control <- if(optim_method == "BFGS"){
      list(
        maxit = maxit,
        factr = reltol,
        fnscale = -1
      )
    }else{
      list(
        maxit = maxit,
        reltol = reltol,
        fnscale = -1
      )
    }
    X <- as.matrix(X)
    if(methods[2] == "RK") stop("For second method, select UNI or PADE (ordering avoided)")
    if(any(dim(X) == 0)) stop("input covariate matrix X, or use fit method instead")
    is_iph <- is(x, "iph")
    EMstep <- eval(parse(text = paste("EMstep_", methods[1], sep = "")))
    if(!all(c(y, rcen) > 0)) stop("data should be positive")
    if(!all(c(weight, rcenweight) >= 0)) stop("weights should be non-negative")
    is_iph <- methods::is(x, "iph")
    if (!is_iph) {
      par_g <- numeric(0)
      inv_g <- function(x) x 
      LL_base <- eval(parse(text = paste("logLikelihoodPH_", methods[2], "s", sep = "")))
      LL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
        ex <- exp(X%*%theta)
        scale1 <- ex[1:length(obs)]
        scale2 <- tail(ex, length(rcens))
        return(LL_base(h, alpha, S, obs, weight, rcens, rcweight, scale1, scale2))
      }
    }else if (is_iph) {
      name <- x@gfun$name
      if(name %in% c("loglogistic", "gev")) stop("not yet available for multi-parameter transforms")
      par_g <- x@gfun$pars
      inv_g <- x@gfun$inverse
      LL_base <- eval(parse(text = paste("logLikelihoodM", x@gfun$name, "_", methods[2], "s", sep = "")))
      LL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
        beta <- theta[1]; B <- theta[2:length(theta)]
        if(beta < 0) return(NA)
        ex <- exp(X%*%B)
        scale1 <- ex[1:length(obs)]
        scale2 <- tail(ex, length(rcens))
        return(LL_base(h, alpha, S, beta, obs, weight, rcens, rcweight, scale1, scale2))
      }
    }

    p <- dim(X)[2]
    n1 <- length(y)
    n2 <- length(rcen)
    ng <- length(par_g)
    
    if(length(weight) == 0) weight <- rep(1, n1)
    if(length(rcenweight) == 0) rcenweight <- rep(1, n2)

    ph_par <- x@pars
    alpha_fit <- clone_vector(ph_par$alpha)
    S_fit <- clone_matrix(ph_par$S)
    if(length(B0) == 0){B_fit <- rep(0, p)
    }else{B_fit <- B0}
    for (k in 1:stepsEM) {
      prop <- exp(X%*%B_fit)
      
      trans <- prop[1:n1] * inv_g(par_g, y)
      trans_cens <- prop[(n1 + 1):(n1 + n2)] * inv_g(par_g, rcen)
      
      A <- data_aggregation(trans, weight)
      B <- data_aggregation(trans_cens, rcenweight)
      
      epsilon1 <- switch(which(methods[1] == c("RK", "UNI","PADE")),
                         if(!is.na(rkstep)){rkstep} else{default_step_length(S_fit)},
                         if(!is.na(uni_epsilon)){uni_epsilon} else{1e-4},
                         0)
      epsilon2 <- switch(which(methods[2] == c("RK", "UNI","PADE")),
                         if(!is.na(rkstep)){rkstep} else{default_step_length(S_fit)},
                         if(!is.na(uni_epsilon)){uni_epsilon} else{1e-4},
                         0)
      EMstep(epsilon1, alpha_fit, S_fit, A$un_obs, A$weights, B$un_obs, B$weights)
      theta <- c(par_g, B_fit)
      opt <- suppressWarnings(optim(par = theta, fn = LL, 
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
                                    )
                              )
      par_g <- head(opt$par, ng)
      B_fit <- tail(opt$par, p)
      if (k %% every == 0) {
        cat("\r", "iteration:", k,
            ", logLik:", opt$value,
            sep = " ")
      }
    }
    cat("\n", sep = "")
    x@pars$alpha <- alpha_fit
    x@pars$S <- S_fit
    x@fit <- list(
      logLik = opt$value,
      nobs = sum(A$weights),
      hessian = opt$hessian
    )
    s <- sph(x, type = "reg")
    s@gfun$pars <- par_g
    s@coefs$B <- B_fit
    return(s)
  }
)


# reg_g_specs <- function(name){
#   if(name == "Homogeneous"){
#     inv_g <- function(t, w, beta) return(list(obs = t, weight = w)) 
#     mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
#       ex <- exp(X%*%theta)
#       scale1 <- ex[1:length(obs)]
#       scale2 <- tail(ex, length(rcens))
#       #o1 <- order(scale1 * obs)
#       #o2 <- order(scale2 * rcens)
#       #return(- logLikelihoodPH_RKs(h, alpha, S, obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
#       return(- logLikelihoodPH_UNIs(0.0001, alpha, S, obs, weight, rcens, rcweight, scale1, scale2))
#       #return(- logLikelihoodPH_PADEs(alpha, S, obs, weight, rcens, rcweight, scale1, scale2))
#     }
#   }
#   else if(name == "weibull"){
#     inv_g <- function(t, w, beta) return(list(obs = t^{beta}, weight = w)) 
#     mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
#       beta <- theta[1]; B <- theta[2:length(theta)]
#       if(beta < 0) return(NA)
#       ex <- exp(X%*%B)
#       scale1 <- ex[1:length(obs)]
#       scale2 <- tail(ex, length(rcens))
#       #o1 <- order(scale1 * inv_g(obs, weight, beta)$obs)
#       #o2 <- order(scale2 * inv_g(rcens, rcweight, beta)$obs)
#       #return(- logLikelihoodMWeib_RKs(h, alpha, S, beta, obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
#       return(- logLikelihoodMweibull_UNIs(0.00001, alpha, S, beta, obs, weight, rcens, rcweight, scale1, scale2))
#       #return(- logLikelihoodMweibull_PADEs(alpha, S, beta, obs, weight, rcens, rcweight, scale1, scale2))
#     }
#   }
#   else if(name == "pareto"){
#     inv_g <- function(t, w, beta) return(list(obs = log(t/beta + 1), weight = w))
#     mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
#       beta <- theta[1]; B <- theta[2:length(theta)]
#       if(beta < 0) return(NA)
#       ex <- exp(X%*%B)
#       scale1 <- ex[1:length(obs)]
#       scale2 <- tail(ex, length(rcens))
#       #o1 <- order(scale1 * inv_g(obs, weight, beta)$obs)
#       #o2 <- order(scale2 * inv_g(rcens, rcweight, beta)$obs)
#       #return(- logLikelihoodMPar_RKs(h, alpha, S, beta, obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
#       return(- logLikelihoodMpareto_UNIs(0.00001, alpha, S, beta, obs, weight, rcens, rcweight, scale1, scale2))
#       #return(- logLikelihoodMpareto_PADEs(alpha, S, beta, obs, weight, rcens, rcweight, scale1, scale2))
#       }
#   }
#   else if(name == "lognormal"){
#     inv_g <- function(t, w, beta) return(list(obs = log(t + 1)^{beta}, weight = w)) 
#     mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
#       beta <- theta[1]; B <- theta[2:length(theta)]
#       if(beta < 0) return(NA)
#       ex <- exp(X%*%B)
#       scale1 <- ex[1:length(obs)]
#       scale2 <- tail(ex, length(rcens))
#       #o1 <- order(scale1 * inv_g(obs, weight, beta)$obs)
#       #o2 <- order(scale2 * inv_g(rcens, rcweight, beta)$obs)
#       #return(- logLikelihoodMLogNormal_RKs(h, alpha, S, beta, obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
#       return(- logLikelihoodMlognormal_UNIs(0.00001, alpha, S, beta, obs, weight, rcens, rcweight, scale1, scale2))
#       #return(- logLikelihoodMlognormal_PADEs(alpha, S, beta, obs, weight, rcens, rcweight, scale1, scale2))
#     }
#   }
#   else if(name == "loglogistic"){
#     inv_g <- function(t, w, beta) return(list(obs = log((t/beta[1])^{beta[2]} + 1), weight = w))
#     mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
#       beta <- theta[1:2]; B <- theta[3:length(theta)]
#       if(beta[1] < 0 | beta[2] < 0) return(NA)
#       ex <- exp(X%*%B)
#       scale1 <- ex[1:length(obs)]
#       scale2 <- tail(ex, length(rcens))
#       #o1 <- order(scale1 * inv_g(obs, weight, beta)$obs)
#       #o2 <- order(scale2 * inv_g(rcens, rcweight, beta)$obs)
#       #return(- logLikelihoodMLogLogistic_RKs(h, alpha, S, beta, obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
#       return(- logLikelihoodMloglogistic_UNIs(0.00001, alpha, S, beta, obs, weight, rcens, rcweight, scale1, scale2))
#       #return(- logLikelihoodMloglogistic_PADEs(alpha, S, beta, obs, weight, rcens, rcweight, scale1, scale2))
#     }
#   }
#   else if(name == "gompertz"){
#     inv_g <- function(t, w, beta) return(list(obs = (exp(t * beta) - 1) / beta, weight = w))
#     mLL <- function(h, alpha, S, theta, obs, weight, rcens, rcweight, X) {
#       beta <- theta[1]; B <- theta[2:length(theta)]
#       if(beta < 0) return(NA)
#       ex <- exp(X%*%B)
#       scale1 <- ex[1:length(obs)]
#       scale2 <- tail(ex, length(rcens))
#       #cat("\r", "Beta:", beta, sep = " ")
#       #o1 <- order(scale1 * inv_g(obs, weight, beta)$obs)
#       #o2 <- order(scale2 * inv_g(rcens, rcweight, beta)$obs)
#       #return(- logLikelihoodMGomp_RKs(h, alpha, S, beta, obs[o1], weight[o1], rcens[o2], rcweight[o2], scale1[o1], scale2[o2]))
#       return(- logLikelihoodMgompertz_UNIs(0.00001, alpha, S, beta, obs, weight, rcens, rcweight, scale1, scale2))
#       #return(- logLikelihoodMgompertz_PADEs(alpha, S, beta, obs, weight, rcens, rcweight, scale1, scale2))
#     }
#   }
#   else{
#     stop("fit for this gfun is not yet implemented")
#   }
#   return(list(inv_g = inv_g, mLL = mLL))
# }
# 

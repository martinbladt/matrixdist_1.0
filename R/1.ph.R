#' Phase Type distributions
#'
#' Class of objects for phase type distributions
#'
#' @slot name name of the phase type distribution
#' @slot pars a list comprising of the parameters.
#'
#' @return
#' @export
#'
setClass("ph",
         slots = list(
           name = "character",
           pars = "list",
           fit = "list"
         ),
         prototype = list(
           name = NA_character_,
           pars = list(),
           fit = list()
         )
)

#' Constructor Function for phase type distributions
#'
#' @param alpha a probability vector.
#' @param S a sub-intensity matrix.
#' @param structure a valid ph structure
#' @param dimension the dimension of the ph structure (if provided)
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
ph <- function(alpha = NULL, S = NULL, structure = NULL, dimension = 3) {
  if (any(is.null(alpha)) & any(is.null(S)) & is.null(structure)) {
    stop("input a vector and matrix, or a structure")
  }
  if (!is.null(structure)) {
    rs <- random_structure(dimension, structure = structure)
    alpha <- rs[[1]]
    S <- rs[[2]]
    name <- structure
  } else {
    if (dim(S)[1] != dim(S)[2]) {
      stop("matrix T should be square")
    }
    if (length(alpha) != dim(S)[1]) {
      stop("incompatible dimensions")
    }
    name <- "custom"
  }
  new("ph",
      name = paste(name, " ph(", length(alpha), ")", sep = ""),
      pars = list(alpha = alpha, S = S)
  )
}

#' Show Method for phase type distributions
#'
#' @param x an object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#'
setMethod("show", "ph", function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
})

#' Simulation Method for phase type distributions
#'
#' @param x an object of class \linkS4class{ph}.
#' @param n length of realization.
#'
#' @return A realization of a phase type data
#' @export
#'
#' @examples
#'
setMethod("sim", c(x = "ph"), function(x, n = 1000) {
  U <- rphasetype(n, x@pars$alpha, x@pars$S)
  return(U)
})

#' Density Method for phase type distributions
#'
#' @param x an object of class \linkS4class{ph}.
#' @param y locations
#'
#' @return Density evaluated at locations
#' @export
#'
#' @examples
#'
setMethod("dens", c(x = "ph"), function(x, y = seq(0, quan(x, .95)$quantile, length.out = 10)) {
  y_inf <- (y == Inf)
  dens <- y
  dens[!y_inf] <- phdensity(y, x@pars$alpha, x@pars$S)
  dens[y_inf] <- 0
  return(list(y = y, dens = dens))
})

#' Distribution Method for phase type distributions
#'
#' @param x an object of class \linkS4class{ph}.
#' @param q locations
#'
#' @return CDF evaluated at locations
#' @export
#'
#' @examples
#'
setMethod("cdf", c(x = "ph"), function(x, 
                                     q = seq(0, quan(x, .95)$quantile, length.out = 10),
                                     lower.tail = TRUE) {
  q_inf <- (q == Inf)
  cdf <- q
  cdf[!q_inf] <- phcdf(q[!q_inf], x@pars$alpha, x@pars$S, lower.tail)
  cdf[q_inf] <- as.numeric(1 * lower.tail)
  return(list(q = q, cdf = cdf))
})

#' Hazard rate Method for phase type distributions
#'
#' @param x an object of class \linkS4class{ph}.
#' @param y locations
#'
#' @return Hazard rate evaluated at locations
#' @export
#'
#' @examples
#'
setMethod("haz", c(x = "ph"), function(x, y = seq(0, quan(x, .95)$quantile, length.out = 10)) {
  d <- dens(x, y)$dens
  s <- cdf(x, y, lower.tail = FALSE)$cdf
  return(list(y = y, haz = d/s))
})

#' Quantile Method for phase type distributions
#'
#' @param x an object of class \linkS4class{ph}.
#' @param p probabilities
#'
#' @return quantiles
#' @export
#'
#' @examples
#'
setMethod("quan", c(x = "ph"), function(x, 
                                     p = seq(0, 1, length.out = 10)) {
  quan <- numeric(length(p))
  for(i in seq_along(p)){
    quan[i] <- uniroot(f = function(q) p[i] - cdf(x, 1/(1 - q) - 1)$cdf, interval = c(0, 1))$root
  }
  return(list(p = p, quantile = 1/(1 - quan) - 1))
})

#' Fit Method for ph Class
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
  "fit", c(x = "ph", y = "ANY"),
  function(x,
           y,
           weight = numeric(0),
           rcen = numeric(0),
           rcenweight = numeric(0),
           stepsEM = 1000) {
    is_iph <- is(x, "iph")
    if(is_iph){
      name <- x@gfun$name
      par_g <- x@gfun$pars
      specs <- g_specs(name) 
      inv_g <- specs$inv_g 
      mLL <- specs$mLL
    }
    A <- data_aggregation(y, weight); y <- A$un_obs; weight <- A$weights
    B <- data_aggregation(rcen, rcenweight); rcen <- B$un_obs; rcenweight <- B$weights
    
    ph_par <- x@pars
    pi_fit <- clone_vector(ph_par$alpha)
    T_fit <- clone_matrix(ph_par$S)
    
    if(!is_iph){
      for (k in 1:stepsEM) {
        RKstep <- default_step_length(T_fit)
        EMstep_RK(RKstep, pi_fit, T_fit, y, weight, rcen, rcenweight)
        if (k %% 100 == 0) {
          cat("\r", "iteration:", k,
              ", logLik:", logLikelihoodPH_RK(RKstep, pi_fit, T_fit, y, weight, rcen, rcenweight),
              sep = " ")
        }
      }
      cat("\n", sep = "")
      x@pars$alpha <- pi_fit
      x@pars$S <- T_fit
    }
    if(is_iph){
      for (k in 1:stepsEM) {
        trans <- inv_g(y, weight, par_g)
        trans_cens <- inv_g(rcen, rcenweight, par_g)
        RKstep <- default_step_length(T_fit)
        EMstep_RK(RKstep, pi_fit, T_fit, trans$obs, trans$weight, trans_cens$obs, trans_cens$weight)
        opt <- suppressWarnings(
          optim(par = par_g,
                fn = mLL,
                h = RKstep,
                alpha = pi_fit,
                S = T_fit, 
                obs = y, 
                weight = weight, 
                rcens = rcen, 
                rcweight = rcenweight,
                hessian = (k == stepsEM),
                method = ifelse(k == stepsEM, "Nelder-Mead", "Nelder-Mead"),
                control = list(
                  maxit = ifelse(k == stepsEM, 1000, 50),
                  reltol = ifelse(k == stepsEM, 1e-8, 1e-6))
          )
          )
        par_g <- opt$par
        if (k %% 10 == 0) {
          cat("\r", "iteration:", k,
              ", logLik:", - opt$value,
              sep = " ")
        }
      }
      cat("\n", sep = "")
      x@pars$alpha <- pi_fit
      x@pars$S <- T_fit
      x@fit <- list(cov = safe_cov(opt$hessian),
                    loglik = - opt$value)
      x <- iph(x, gfun = name, gfun_pars = par_g)
    }
    return(x)
  }
)

data_aggregation <- function(y, w){
  y <- sort(as.numeric(y))
  un_obs <- unique(y)
  if (length(w) == 0) {
    w <- rep(1, length(y))
  }
  observations <- cbind(y, w)
  mat <- data.frame(observations)
  names(mat) <- c("obs", "weight")
  cum_weight <- numeric(0)
  for (i in un_obs) {
    cum_weight <- c(cum_weight, sum(mat$weight[which(mat$obs == i)]))
  }
  return(list(un_obs = un_obs, weights = cum_weight))
}


g_specs <- function(name){
  if(name == "Weibull"){
    inv_g <- function(t, w, beta) return(list(obs = t^{beta}, weight = w)) 
    mLL <- function(h, alpha, S, beta, obs, weight, rcens, rcweight) {
      if(beta < 0) return(NA)
      return(- logLikelihoodMWeib_RK(h, alpha, S, beta, obs, weight, rcens, rcweight))
    }
  }
  else if(name == "Pareto"){
    inv_g <- function(t, w, beta) return(list(obs = log(t/beta + 1), weight = w))
    mLL <- function(h, alpha, S, beta, obs, weight, rcens, rcweight) {
      if(beta < 0) return(NA)
      return(- logLikelihoodMPar_RK(h, alpha, S, beta, obs, weight, rcens, rcweight))
    }
  }
  else if(name == "LogNormal"){
    inv_g <- function(t, w, beta) return(list(obs = log(t/beta + 1), weight = w))
    mLL <- function(h, alpha, S, beta, obs, weight, rcens, rcweight) {
      if(beta < 0) return(NA)
      return(- logLikelihoodMLogNormal_RK(h, alpha, S, beta, obs, weight, rcens, rcweight))
    }
  }
  else if(name == "LogLogistic"){
    inv_g <- function(t, w, beta) return(list(obs = log((t/beta[1])^{beta[2]} + 1), weight = w))
    mLL <- function(h, alpha, S, beta, obs, weight, rcens, rcweight) {
      if(beta[1] < 0 | beta[2] < 0) return(NA)
      return(- logLikelihoodMLogLogistic_RK(h, alpha, S, beta, obs, weight, rcens, rcweight))
    }
  }
  else if(name == "Gompertz"){
    inv_g <- function(t, w, beta) return(list(obs = (exp(t * beta) - 1) / beta, weight = w))
    mLL <- function(h, alpha, S, beta, obs, weight, rcens, rcweight) {
      if(beta < 0) return(NA)
      return(- logLikelihoodMGomp_RK(h, alpha, S, beta, obs, weight, rcens, rcweight))
    }
  }
  else if(name == "GEVD"){
    inv_g <- reversTransformData
    mLL <- function(h, alpha, S, beta, obs, weight, rcens, rcweight) {
      if(beta[2] < 0) return(NA)
      return(- logLikelihoodMGEV_RK(h, alpha, S, beta, obs, weight, rcens, rcweight))
    }
  }else{
    stop("fit for this gfun is not yet implemented")
  }
  return(list(inv_g = inv_g, mLL = mLL))
}

#' Calculate Covariance Matrix Safely
#'
#' @param hess a Hessian matrix from a model fit.
#'
#' @return Fisher information (estimated covariance matrix)
#' @export
#'
safe_cov <- function(hess) {
  hessinverse <- tryCatch(solve(hess), error = function(e) {
    warning("Hessian can't be inverted")
    return(matrix(NA, nrow = nrow(hess), ncol = ncol(hess)))
  })
  hessinverse
}


#' Coef Method for ph Class
#'
#' @param object an object of class \linkS4class{ph}.
#'
#' @return parameters of ph model
#' @export
#'
setMethod("coef", c(object = "ph"), function(object) {
  object@pars
})

#' Plot Method for phase type distributions
#'
#' @param x an object of class \linkS4class{ph}.
#' @param y a dataset
#'
#' @export
#'
#' @examples
#'
setMethod("plot", c(x = "ph"), function(x, y = NULL) {
  if (all(is.null(y))) {
    sq <- seq(1e-20, quan(x, 0.99)$quan, length.out = 1000)
    phd <- phdensity(sq, x@pars$alpha, x@pars$S)
    base::plot(sq, phd, type = "l", xlab = "y", ylab = "density")
  }
  if (!all(is.null(y))) {
    sq <- seq(1e-20, max(y), length.out = 1000)
    phd <- phdensity(sq, x@pars$alpha, x@pars$S)
    mx_h <- max(hist(y, breaks = 100, plot = FALSE)$density)
    mx_d <- max(phd)
    hist(y, breaks = 100, freq = FALSE, ylim = c(0, max(mx_h, mx_d)))
    lines(sq, phd, col = "red")
  }
})

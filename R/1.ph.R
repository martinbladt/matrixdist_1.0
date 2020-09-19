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
           pars = "list"
         ),
         prototype = list(
           name = NA_character_,
           pars = list()
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
    name <- "Custom"
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
  if (length(object@name) > 0) {
    cat("phase-type name: ", object@name, "\n", sep = "")
    cat("parameters: ", "\n", sep = "")
    print(object@pars)
  } else {
    return()
  }
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
setMethod("r", c(x = "ph"), function(x, n = 1000) {
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
setMethod("d", c(x = "ph"), function(x, y = seq(0, 5, length.out = 100)) {
  dens <- phdensity(y, x@pars$alpha, x@pars$S)
  return(cbind(y = y, dens = dens))
})

#' Distribution Method for phase type distributions
#'
#' @param x an object of class \linkS4class{ph}.
#' @param y locations
#'
#' @return CDF evaluated at locations
#' @export
#'
#' @examples
#'
setMethod("p", c(x = "ph"), function(x, 
                                     q = seq(0, 5, length.out = 100),
                                     lower.tail = TRUE) {
  cdf <- phcdf(q, x@pars$alpha, x@pars$S, lower.tail)
  return(cbind(q = q, cdf = cdf))
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
           xweight = numeric(0),
           rcen = numeric(0),
           rcenweight = numeric(0),
           stepsEM = 1000) {
    is_iph <- is(x, "iph")
    if(is_iph){
      name <- x@gfun$name
      par_g <- x@gfun$pars
      x <- x@ph
      if(name == "Weibull"){
        inv_g <- function(t, beta) t^{beta}
        mLL <- function(h, alpha, S, beta, obs, weight, rcens, rcweight) {
          if(beta < 0) return(NA)
          return(- logLikelihoodMWeib_RK(h, alpha, S, beta, obs, weight, rcens, rcweight))
        }
      }
      else if(name == "Pareto"){
        inv_g <- function(t, beta) log(t/beta + 1)
        mLL <- function(h, alpha, S, beta, obs, weight, rcens, rcweight) {
          if(beta < 0) return(NA)
          return(- logLikelihoodMPar_RK(h, alpha, S, beta, obs, weight, rcens, rcweight))
        }
      }
      else if(name == "Gompertz"){
        inv_g <- function(t, beta) (exp(t * beta) - 1) / beta
        mLL <- function(h, alpha, S, beta, obs, weight, rcens, rcweight) {
          if(beta < 0) return(NA)
          return(- logLikelihoodMGomp_RK(h, alpha, S, beta, obs, weight, rcens, rcweight))
        }
      }
      else if(name == "GEVD"){
        inv_g <- function(t, beta){
          mu <- beta[1]; sigma <- beta[2]; xi <- beta[3]
          if(xi == 0) return(exp(-(t - mu) / sigma))
          else return((1 + (xi / sigma) * (t - mu))^{-1 / xi})
        }
        mLL <- function(h, alpha, S, beta, obs, weight, rcens, rcweight) {
          if(beta[2] < 0) return(NA)
          return(- logLikelihoodMGEV_RK(h, alpha, S, beta, obs, weight, rcens, rcweight))
        }
      }else{
        stop("fit for this gfun is not yet implemented")
      }
    }
    y <- sort(as.numeric(y))
    un_obs <- unique(y)
    if (min(y) <= 0) {
      stop("data should be positive")
    }
    if (length(xweight) == 0) {
      xweight <- rep(1, length(y))
    }
    observations <- cbind(y, xweight)
    mat <- data.frame(observations)
    names(mat) <- c("obs", "weight")
    cum_weight <- NULL
    for (i in un_obs) {
      cum_weight <- c(cum_weight, sum(mat$weight[which(mat$obs == i)]))
    }
    ph_par <- x@pars
    pi_fit <- clone_vector(ph_par$alpha)
    T_fit <- clone_matrix(ph_par$S)
    
    if(!is_iph){
      for (k in 1:stepsEM) {
        RKstep <- default_step_length(T_fit)
        EMstep_RK(RKstep, pi_fit, T_fit, un_obs, cum_weight, rcen, rcenweight)
        if (k %% 100 == 0) {
          cat("\r", "iteration:", k,
              ", logLik:", logLikelihoodPH_RK(RKstep, pi_fit, T_fit, un_obs, cum_weight, rcen, rcenweight),
              sep = " ")
        }
      }
      cat("\n", sep = "")
      x@pars$alpha <- pi_fit
      x@pars$S <- T_fit
    }
    if(is_iph){
      for (k in 1:stepsEM) {
        trans_obs <- inv_g(un_obs, par_g)
        RKstep <- default_step_length(T_fit)
        EMstep_RK(RKstep, pi_fit, T_fit, trans_obs, cum_weight, rcen, rcenweight)
        opt <- suppressWarnings(optim(par = par_g, fn = mLL, h = RKstep, alpha = pi_fit, S = T_fit, obs = un_obs, weight = cum_weight, rcens = rcen, rcweight = rcenweight))
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
      x <- iph(x, gfun = name, gfun_pars = par_g)
    }
    return(x)
  }
)

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
setMethod("m_plot", c(x = "ph"), function(x, y = NULL) {
  if (all(is.null(y))) {
    sq <- seq(1e-20, 5, length.out = 1000)
    phd <- phdensity(sq, x@pars$alpha, x@pars$S)
    plot(sq, phd, type = "l", xlab = "y", ylab = "density")
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

#' Phase Type distributions
#'
#' Class of objects for phase type distributions
#'
#' @slot name name of the phase type distribution
#' @slot pars a list comprising of the parameters
#' @slot fit a list containing estimation information
#'
#' @return Class object
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
#' @param dimension the dimension of the ph structure (if structure is provided)
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' ph(structure = "gcoxian", dim = 5)
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
      stop("matrix S should be square")
    }
    if (length(alpha) != dim(S)[1]) {
      stop("incompatible dimensions")
    }
    name <- "custom"
  }
  methods::new("ph",
    name = paste(name, " ph(", length(alpha), ")", sep = ""),
    pars = list(alpha = alpha, S = S)
  )
}

#' Sum Method for phase type distributions
#'
#' @param e1 an object of class \linkS4class{ph}.
#' @param e2 an object of class \linkS4class{ph}.
#' @export
#'
setMethod("+", signature(e1 = "ph", e2 = "ph"), 
          function (e1, e2){
            if(methods::is(e1, "iph") | methods::is(e2, "iph")) stop("objects to be added should be ph")
            L <- sumPH(e1@pars$alpha, e1@pars$S, e2@pars$alpha, e2@pars$S)
            return(ph(alpha = L$pi, S = L$T))
          }
          )

kronecker_sum <- function(A, B){
  n <- nrow(A); m <- nrow(B)
  kronecker(A, diag(m)) + kronecker(diag(n), B)
}

#' Minimum Method for phase type distributions
#'
#' @param x1 an object of class \linkS4class{ph}.
#' @param x2 an object of class \linkS4class{ph}.
#' 
#' @return An object of class \linkS4class{ph}.
#' @export
#'
setMethod("minimum", signature(x1 = "ph", x2 = "ph"), 
          function (x1, x2){
            alpha <- kronecker(x1@pars$alpha, x2@pars$alpha)
            S <- kronecker_sum(x1@pars$S, x2@pars$S)
            return(ph(alpha = alpha, S = S))
          }
)

#' Maximum Method for phase type distributions
#'
#' @param x1 an object of class \linkS4class{ph}.
#' @param x2 an object of class \linkS4class{ph}.
#' @export
#'
setMethod("maximum", signature(x1 = "ph", x2 = "ph"), 
          function (x1, x2){
            n1 <- length(x1@pars$alpha)
            n2 <- length(x2@pars$alpha)
            alpha <- c(kronecker(x1@pars$alpha, x2@pars$alpha), rep(0, n1 + n2))
            S1 <- rbind(kronecker_sum(x1@pars$S, x2@pars$S), matrix(0, n1 + n2, n1 * n2))
            S2 <- rbind(kronecker(diag(n1), -rowSums(x2@pars$S)), x1@pars$S, matrix(0, n2, n1))
            S3 <- rbind(kronecker(-rowSums(x1@pars$S), diag(n2)), matrix(0, n1, n2), x2@pars$S)
            return(ph(alpha = alpha, S = cbind(S1, S2, S3)))
          }
)

#' Moment Method for phase type distributions
#'
#' @param x an object of class \linkS4class{ph}.
#' @param k a positive integer.
#' @export
#'
setMethod("moment", signature(x = "ph"), 
          function (x, k = 1){
            if(k%%1 != 0 | k <= 0) return("k should be a positive integer")
            if(methods::is(x, "iph")) warning("moment of undelying ph structure is provided for iph objects")
            return(phmoment(k, x@pars$alpha, x@pars$S))
          }
)


#' Show Method for phase type distributions
#'
#' @param object an object of class \linkS4class{ph}.
#' @export
#'
setMethod("show", "ph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
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
#' obj <- ph(structure = "general")
#' sim(obj, n = 100)
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
#' obj <- ph(structure = "general")
#' dens(obj, c(1, 2, 3))
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
#' @param lower.tail cdf(TRUE) or tail(FALSE)
#'
#' @return CDF evaluated at locations
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' cdf(obj, c(1, 2, 3))
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
#' obj <- ph(structure = "general")
#' haz(obj, c(1, 2, 3))
setMethod("haz", c(x = "ph"), function(x, y = seq(0, quan(x, .95)$quantile, length.out = 10)) {
  d <- dens(x, y)$dens
  s <- cdf(x, y, lower.tail = FALSE)$cdf
  return(list(y = y, haz = d / s))
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
#' obj <- ph(structure = "general")
#' quan(obj, c(0.5, 0.9, 0.99))
setMethod("quan", c(x = "ph"), function(x,
                                        p = seq(0, 1, length.out = 10)) {
  quan <- numeric(length(p))
  for (i in seq_along(p)) {
    quan[i] <- stats::uniroot(f = function(q) p[i] - cdf(x, 1 / (1 - q) - 1)$cdf, interval = c(0, 1))$root
  }
  return(list(p = p, quantile = 1 / (1 - quan) - 1))
})

#' Fit Method for ph Class
#'
#' @param x an object of class \linkS4class{ph}.
#' @param y vector or data.
#' @param weight vector of weights.
#' @param rcen vector of right-censored observations
#' @param rcenweight vector of weights for right-censored observations.
#' @param stepsEM number of EM steps to be performed.
#' @param rkstep Runge-Kutta step size (optional)
#' @param maxit maximum number of iterations when optimizing g function.
#' @param reltol relative tolerance when optimizing g function.
#' @param every number of iterations between likelihood display updates.
#' 
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' data <- sim(obj)
#' fit(obj, data)
setMethod(
  "fit", c(x = "ph", y = "ANY"),
  function(x,
           y,
           weight = numeric(0),
           rcen = numeric(0),
           rcenweight = numeric(0),
           stepsEM = 1000,
           rkstep = NA,
           maxit = 100,
           reltol = 1e-8,
           every = 100) {
    is_iph <- methods::is(x, "iph")
    if (is_iph) {
      par_g <- x@gfun$pars
      inv_g <- x@gfun$inverse
      mLL <- eval(parse(text = paste("logLikelihoodM", x@gfun$name, "_RK", sep = "")))
    }
    A <- data_aggregation(y, weight)
    y <- A$un_obs
    weight <- A$weights
    B <- data_aggregation(rcen, rcenweight)
    rcen <- B$un_obs
    rcenweight <- B$weights

    ph_par <- x@pars
    pi_fit <- clone_vector(ph_par$alpha)
    T_fit <- clone_matrix(ph_par$S)

    if (!is_iph) {
      for (k in 1:stepsEM) {
        if(!is.na(rkstep)) RKstep <- rkstep else  RKstep <- default_step_length(T_fit)
        EMstep_RK(RKstep, pi_fit, T_fit, y, weight, rcen, rcenweight)
        if (k %% every == 0) {
          cat("\r", "iteration:", k,
            ", logLik:", logLikelihoodPH_RK(RKstep, pi_fit, T_fit, y, weight, rcen, rcenweight),
            sep = " "
          )
        }
      }
      cat("\n", sep = "")
      x@pars$alpha <- pi_fit
      x@pars$S <- T_fit
    }
    if (is_iph) {
      trans_weight <- weight 
      trans_rcenweight <- rcenweight
      for (k in 1:stepsEM) {
        if(x@gfun$name != "gev") {trans <- inv_g(par_g, y); trans_cens <- inv_g(par_g, rcen)
        }else{ t <- inv_g(par_g, y, weight); tc <- inv_g(par_g, rcen, rcenweight) 
        trans <- t$obs; trans_weight <- t$weight; trans_cens <- tc$obs; trans_rcenweight <- tc$weight}
        if(!is.na(rkstep)) RKstep <- rkstep else  RKstep <- default_step_length(T_fit)
        EMstep_RK(RKstep, pi_fit, T_fit, trans, trans_weight, trans_cens, trans_rcenweight)
        opt <- suppressWarnings(
          stats::optim(
            par = par_g,
            fn = mLL,
            h = RKstep,
            pi = pi_fit,
            T = T_fit,
            obs = y,
            weight = weight,
            rcens = rcen,
            rcweight = rcenweight,
            hessian = (k == stepsEM),
            control = list(
              maxit = maxit,
              reltol = reltol
            )
          )
        )
        par_g <- opt$par
        if (k %% every == 0) {
          cat("\r", "iteration:", k,
            ", logLik:", -opt$value,
            sep = " "
          )
        }
      }
      cat("\n", sep = "")
      x@pars$alpha <- pi_fit
      x@pars$S <- T_fit
      x@fit <- list(
        cov = safe_cov(opt$hessian),
        loglik = -opt$value
      )
      x <- iph(x, gfun = x@gfun$name, gfun_pars = par_g)
    }
    return(x)
  }
)

data_aggregation <- function(y, w) {
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

#' Calculate Covariance Matrix Safely
#'
#' @param hess a Hessian matrix from a model fit.
#'
#' @return Fisher information (estimated covariance matrix)
#'
safe_cov <- function(hess) {
  hessinverse <- tryCatch(solve(hess), error = function(e) {
    warning("hessian can't be inverted")
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
#' @examples
#' obj <- ph(structure = "general")
#' coef(obj)
setMethod("coef", c(object = "ph"), function(object) {
  object@pars
})

#' Inhomogeneous Phase Type distributions
#'
#' Class of objects for inhomogeneous phase type distributions
#'
#' @slot name name of the phase type distribution
#' @slot pars a list comprising of the parameters.
#'
#' @return
#' @export
#'
setClass("iph",
  contains = c("ph"),
  slots = list(
    gfun = "list",
    scale = "numeric"
  )
)

#' Constructor Function for inhomogeneous phase type distributions
#'
#' @param alpha a probability vector.
#' @param S a sub-intensity matrix.
#' @param structure a valid ph structure
#' @param dimension the dimension of the ph structure (if provided)
#' @param fun inhomogeneity transform
#' @param funpars the parameters of the inhomogeneity function
#'
#' @return An object of class \linkS4class{iph}.
#' @export
#'
#' @examples
iph <- function(ph = NULL, gfun = NULL, gfun_pars = NULL, alpha = NULL, S = NULL, structure = NULL, dimension = 3) {
  if (all(is.null(c(gfun, gfun_pars)))) {
    stop("input inhomogeneity function and parameters")
  }
  if (is.null(ph)) {
    ph <- ph(alpha = alpha, S = S, structure = structure, dimension = dimension)
  }
  if (!gfun %in% c("Pareto", "Weibull", "LogLogistic", "Gompertz", "GEVD")) {
    stop("invalid gfun")
  }
  if (gfun %in% c("Pareto", "Weibull", "Gompertz")) {
    if (length(gfun_pars) != 1 | sum(gfun_pars <= 0) > 0) {
      stop("gfun parameter should be positive and of length one")
    } else {
      names(gfun_pars) <- "beta"
    }
  }
  if (gfun %in% c("GEVD")) {
    if (length(gfun_pars) != 3 | (gfun_pars[2] > 0) == FALSE) {
      stop("gfun parameter should be of length three: mu, sigma, xi, and sigma > 0")
    } else {
      names(gfun_pars) <- c("mu", "sigma", "xi")
    }
  }
  if (gfun %in% c("LogLogistic")) {
    if (length(gfun_pars) != 2 | (gfun_pars[1] <= 0) | (gfun_pars[2] <= 0)) {
      stop("gfun parameter should be positive and of length two: alpha, theta > 0")
    } else {
      names(gfun_pars) <- c("alpha", "theta")
    }
  }
  new("iph",
    name = paste("inhomogeneous ", ph@name, sep = ""),
    pars = ph@pars,
    gfun = list(name = gfun, pars = gfun_pars),
    scale = 1
  )
}

#' Show Method for inhomogeneous phase type distributions
#'
#' @param x an object of class \linkS4class{iph}.
#' @export
#'
#' @examples
#'
setMethod("show", "iph", function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("g-function name: ", object@gfun$name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@gfun$pars)
  cat("scale: ", "\n", sep = "")
  print(object@scale)
})

#' Simulation Method for inhomogeneous phase type distributions
#'
#' @param x an object of class \linkS4class{iph}.
#' @param n length of realization.
#'
#' @return A realization of inhomogeneous phase type data
#' @export
#'
#' @examples
#'
setMethod("sim", c(x = "iph"), function(x, n = 1000) {
  name <- x@gfun$name
  pars <- x@gfun$pars
  scale <- x@scale
  if (name %in% c("Pareto", "Weibull", "LogLogistic","Gompertz")) {
    U <- scale * riph(n, name, x@pars$alpha, x@pars$S, pars)
  }
  if (name %in% c("GEVD")) {
    U <- scale * rmatrixGEVD(n, x@pars$alpha, x@pars$S, pars[1], pars[2], pars[3])
  }
  return(U)
})

#' Density Method for inhomogeneous phase type distributions
#'
#' @param x an object of class \linkS4class{iph}.
#' @param y locations
#'
#' @return Density evaluated at locations
#' @export
#'
#' @examples
#'
setMethod("dens", c(x = "iph"), function(x, y = seq(0, quan(x, .95)$quantile, length.out = 10)) {
  fn <- eval(parse(text = paste("m", x@gfun$name, "den", sep = "")))
  scale <- x@scale
  y_inf <- (y == Inf)
  dens <- y
  dens[!y_inf] <- fn(y/scale, x@pars$alpha, x@pars$S, x@gfun$pars)/scale
  dens[y_inf] <- 0
  return(list(y = y, dens = dens))
})

#' Distribution Method for inhomogeneous phase type distributions
#'
#' @param x an object of class \linkS4class{iph}.
#' @param y locations
#'
#' @return CDF evaluated at locations
#' @export
#'
#' @examples
#'
setMethod("cdf", c(x = "iph"), function(x, 
                                      q = seq(0, quan(x, .95)$quantile, length.out = 10),
                                      lower.tail = TRUE) {
  fn <- eval(parse(text = paste("m", x@gfun$name, "cdf", sep = "")))
  scale <- x@scale
  q_inf <- (q == Inf)
  cdf <- q
  cdf[!q_inf] <- fn(q/scale, x@pars$alpha, x@pars$S, x@gfun$pars, lower.tail)
  cdf[q_inf] <- as.numeric(1 * lower.tail)
  return(list(q = q, cdf = cdf))
})


#' Coef Method for iph Class
#'
#' @param object an object of class \linkS4class{iph}.
#'
#' @return parameters of ph model
#' @export
#'
setMethod("coef", c(object = "iph"), function(object) {
  L <- append(object@pars, unname(object@gfun$pars))
  names(L)[3] <- names(object@gfun$pars)
  L
})


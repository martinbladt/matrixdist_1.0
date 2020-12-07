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
#' @param gfun inhomogeneity transform
#' @param gfunpars the parameters of the inhomogeneity function
#'
#' @return An object of class \linkS4class{iph}.
#' @export
#'
#' @examples
iph <- function(ph = NULL, gfun = NULL, gfun_pars = NULL, alpha = NULL, S = NULL, structure = NULL, dimension = 3, scale = 1) {
  if (all(is.null(c(gfun, gfun_pars)))) {
    stop("input inhomogeneity function and parameters")
  }
  if (is.null(ph)) {
    ph <- ph(alpha = alpha, S = S, structure = structure, dimension = dimension)
  }
  if (!gfun %in% c("Pareto", "Weibull", "LogNormal", "LogLogistic", "Gompertz", "GEVD", "Identity")) {
    stop("invalid gfun")
  }
  if (gfun %in% c("Pareto", "Weibull", "LogNormal", "Gompertz")) {
    if(is.null(gfun_pars))gfun_pars <- 1
    if (length(gfun_pars) != 1 | sum(gfun_pars <= 0) > 0) {
      stop("gfun parameter should be positive and of length one")
    } else {
      names(gfun_pars) <- "beta"
    }
  }
  if (gfun %in% c("GEVD")) {
    if(is.null(gfun_pars))gfun_pars <- c(0, 1, 1)
    if (length(gfun_pars) != 3 | (gfun_pars[2] > 0) == FALSE) {
      stop("gfun parameter should be of length three: mu, sigma, xi, and sigma > 0")
    } else {
      names(gfun_pars) <- c("mu", "sigma", "xi")
    }
  }
  if (gfun %in% c("LogLogistic")) {
    if(is.null(gfun_pars))gfun_pars <- c(1, 1)
    if (length(gfun_pars) != 2 | (gfun_pars[1] <= 0) | (gfun_pars[2] <= 0)) {
      stop("gfun parameter should be positive and of length two: alpha, theta > 0")
    } else {
      names(gfun_pars) <- c("alpha", "theta")
    }
  }
  if (gfun == "Weibull") {
    ginv <- function(beta, t) t^{beta}
  }
  else if (gfun == "Pareto") {
    ginv <- function(beta, t) log(t / beta + 1)
  }
  else if (gfun == "LogNormal") {
    ginv <- function(beta, t) log(t + 1)^{beta}
  }
  else if (gfun == "LogLogistic") {
    ginv <- function(beta, t) log((t / beta[1])^{beta[2]} + 1)
  }
  else if (gfun == "Gompertz") {
    ginv <- function(beta, t) (exp(t * beta) - 1) / beta
  }
  else if (gfun == "GEVD") {
    ginv <- function(beta, t, w) reversTransformData(t, w, beta)
  }
  new("iph",
    name = paste("inhomogeneous ", ph@name, sep = ""),
    pars = ph@pars,
    gfun = list(name = gfun, pars = gfun_pars, inverse = ginv),
    scale = scale,
    fit = ph@fit
  )
}

#' Minimum Method for inhomogeneous phase type distributions
#'
#' @param x1 an object of class \linkS4class{iph}.
#' @param x2 an object of class \linkS4class{iph}.
#' @export
#'
setMethod("minimum", signature(x1 = "iph", x2 = "iph"), 
          function (x1, x2){
            if(x1@gfun$name != x2@gfun$name | !all(x1@gfun$pars == x2@gfun$pars)) stop("g functions (or parameters) are different")
            ph1 <- ph(alpha = x1@pars$alpha, S = x1@pars$S)
            ph2 <- ph(alpha = x2@pars$alpha, S = x2@pars$S)
            return(iph(ph = minimum(ph1, ph2), gfun = x1@gfun$name, gfun_pars = x1@gfun$pars))
          }
)

#' Maximum Method for inhomogeneous phase type distributions
#'
#' @param x1 an object of class \linkS4class{iph}.
#' @param x2 an object of class \linkS4class{iph}.
#' @export
#'
setMethod("maximum", signature(x1 = "iph", x2 = "iph"), 
          function (x1, x2){
            if(x1@gfun$name != x2@gfun$name | !all(x1@gfun$pars == x2@gfun$pars)) stop("g functions (or parameters) are different")
            ph1 <- ph(alpha = x1@pars$alpha, S = x1@pars$S)
            ph2 <- ph(alpha = x2@pars$alpha, S = x2@pars$S)
            return(iph(ph = maximum(ph1, ph2), gfun = x1@gfun$name, gfun_pars = x1@gfun$pars))
          }
)

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
  #cat("scale: ", "\n", sep = "")
  #print(object@scale)
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
  if (name %in% c("Pareto", "Weibull", "LogNormal", "LogLogistic", "Gompertz")) {
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
  fn <- base::eval(parse(text = paste("m", x@gfun$name, "den", sep = "")))
  scale <- x@scale
  y_inf <- (y == Inf)
  dens <- y
  dens[!y_inf] <- fn(y / scale, x@pars$alpha, x@pars$S, x@gfun$pars) / scale
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
  fn <- base::eval(parse(text = paste("m", x@gfun$name, "cdf", sep = "")))
  scale <- x@scale
  q_inf <- (q == Inf)
  cdf <- q
  cdf[!q_inf] <- fn(q / scale, x@pars$alpha, x@pars$S, x@gfun$pars, lower.tail)
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

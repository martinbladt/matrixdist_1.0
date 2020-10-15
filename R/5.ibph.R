#' Inhomogeneous Bivariate Phase Type distributions
#'
#' Class of objects for inhomogeneous bivariate phase type distributions
#'
#' @slot name name of the phase type distribution
#' @slot pars a list comprising of the parameters.
#'
#' @return
#' @export
#'
setClass("ibph",
         contains = c("bph"),
         slots = list(
           gfun = "list"
         )
)

#' Constructor Function for inhomogeneous bivariate phase type distributions
#'
#' @param alpha a probability vector.
#' @param S a sub-intensity matrix.
#' @param structure a valid ph structure
#' @param dimension the dimension of the ph structure (if provided)
#' @param fun inhomogeneity transform
#' @param funpars the parameters of the inhomogeneity function
#'
#' @return An object of class \linkS4class{ibph}.
#' @export
#'
#' @examples
ibph <- function(bph = NULL, gfun = NULL, gfun_pars = NULL, alpha = NULL, T11 = NULL, T12 = NULL, T22 = NULL, structure = NULL, dimensions = c(3, 3)) {
  if (all(is.null(c(gfun, gfun_pars)))) {
    stop("input inhomogeneity function and parameters")
  }
  if (is.null(mph)) {
    mph <- bph(alpha = alpha, T11 = T11, T12 = T12, T22 = T22, structure = structure, dimensions = dimensions)
  }
  if (!gfun %in% c("Pareto", "Weibull")) {
    stop("invalid gfun")
  }
  if (gfun %in% c("Pareto", "Weibull")) {
    if (length(gfun_pars) != 2 | sum(gfun_pars <= 0) > 0) {
      stop("gfun parameter should be positive and of length two")
    } else {
      names(gfun_pars) <- c("beta1", "beta2")
    }
  }
  new("ibph",
      name = paste("inhomogeneous ", bph@name, sep = ""),
      blocks = bph@blocks,
      rew = bph@rew,
      pars = bph@pars,
      gfun = list(name = gfun, pars = gfun_pars)
  )
}

#' Show Method for inhomogeneous bivariate phase type distributions
#'
#' @param x an object of class \linkS4class{ibph}.
#' @export
#'
#' @examples
#'
setMethod("show", "ibph", function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("number of variables: ", dim(object@rew$R)[2], "\n", sep = "")
  cat("reward matrix: ", "\n", sep = "")
  print(object@rew$R)
  cat("g-function name: ", object@gfun$name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@gfun$pars)
})

#' Simulation Method for inhomogeneous bivariate phase type distributions
#'
#' @param x an object of class \linkS4class{ibph}.
#' @param n length of realization.
#'
#' @return A realization of inhomogeneous phase type data
#' @export
#'
#' @examples
#'
setMethod("sim", c(x = "ibph"), function(x, n = 1000) {
  name <- x@gfun$name
  pars <- x@gfun$pars
  if (name %in% c("Pareto", "Weibull")) {
    U <- rimph(n, name, x@pars$alpha, x@pars$S, x@rew$R, pars)
  }
  return(U)
})

#' Density Method for inhomogeneous bivariate phase type distributions
#'
#' @param x an object of class \linkS4class{ibph}.
#' @param y locations
#'
#' @return Density evaluated at locations
#' @export
#'
#' @examples
#'
setMethod("dens", c(x = "ibph"), function(x, y = matrix(c(1, 2, 1, 2), 2, 2)) {
  fn <- base::eval(parse(text = paste("bivm", x@gfun$name, "den", sep = "")))
  dens <- fn(y, x@blocks$alpha0, x@blocks$T11, x@blocks$T12, x@blocks$T22, x@gfun$pars)
  return(list(y = y, dens = dens))
})

#' Distribution Method for inhomogeneous bivariate phase type distributions
#'
#' @param x an object of class \linkS4class{ibph}.
#' @param y locations
#'
#' @return CDF evaluated at locations
#' @export
#'
#' @examples
#'
setMethod("cdf", c(x = "ibph"), function(x, 
                                        q = seq(0, quan(x, .95)$quantile, length.out = 10),
                                        lower.tail = TRUE) {
  #TO BE CONSTRUCTED..
  return(list(q = q, cdf = cdf))
})

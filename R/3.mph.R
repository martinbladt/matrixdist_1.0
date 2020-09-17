#' Multivariate Phase Type distributions
#'
#' Class of objects for multivariate phase type distributions
#'
#' @slot name name of the phase type distribution
#' @slot pars a list comprising of the parameters.
#'
#' @return
#' @export
#'
setClass("mph",
  contains = c("ph"),
  slots = list(
    ph = "ph",
    rew = "list"
  ),
  prototype = list(
    name = NA_character_,
    pars = list()
  )
)

#' Constructor Function for multivariate phase type distributions
#'
#' @param alpha a probability vector.
#' @param S a sub-intensity matrix.
#' @param R a compatible (non-negative) reward matrix
#' @param structure a valid ph structure
#' @param dimension the dimension of the ph structure (if provided)
#'
#' @return An object of class \linkS4class{mph}.
#' @export
#'
#' @examples
mph <- function(ph = NULL, R = NULL, alpha = NULL, S = NULL, structure = NULL, dimension = 3) {
  if (is.null(ph)) {
    ph <- ph(alpha = alpha, S = S, structure = structure, dimension = dimension)
  }
  if (dim(R)[1] != dim(ph@pars$S)[1]) {
    stop("matrix R has wrong dimension: number of rows does not match ph dimension")
  }
  new("mph",
    ph = ph,
    rew = list(R = R)
  )
}

#' Show Method for multivariate phase type distributions
#'
#' @param x an object of class \linkS4class{mph}.
#' @export
#'
#' @examples
#'
setMethod("show", "mph", function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  if (length(object@name) > 0) {
    cat("phase-type name: ", object@ph@name, "\n", sep = "")
    cat("parameters: ", "\n", sep = "")
    print(object@ph@pars)
    cat("number of variables: ", dim(object@rew$R)[2], "\n", sep = "")
    cat("reward matrix: ", "\n", sep = "")
    print(object@rew$R)
  } else {
    return()
  }
})

#' Simulation Method for multivariate phase type distributions
#'
#' @param x an object of class \linkS4class{mph}.
#' @param n length of realization.
#'
#' @return A realization of a phase type data
#' @export
#'
#' @examples
#'
setMethod("r", c(x = "mph"), function(x, n = 1000) {
  U <- rmph(n, x@ph@pars$alpha, x@ph@pars$S, x@rew$R)
  return(U)
})

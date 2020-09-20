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
    rew = "list"
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
mph <- function(ph = NULL, R = NULL, alpha = NULL, S = NULL, structure = NULL, p = 3, dimension = 3, variables = 2) {
  if (is.null(ph)) {
    ph <- ph(alpha = alpha, S = S, structure = structure, dimension = dimension)
  }
  if(all(is.null(R))){R <- random_reward(dim(ph@pars$S)[1], variables)
  }else{
    if (dim(R)[1] != dim(ph@pars$S)[1]) {
      stop("matrix R has wrong dimension: number of rows does not match ph dimension")
  }
  }
  new("mph",
    name = paste("mph based on a ", ph@name, sep = ""),
    pars = ph@pars,
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
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("number of variables: ", dim(object@rew$R)[2], "\n", sep = "")
  cat("reward matrix: ", "\n", sep = "")
  print(object@rew$R)
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
setMethod("sim", c(x = "mph"), function(x, n = 1000) {
  U <- rmph(n, x@pars$alpha, x@pars$S, x@rew$R)
  return(U)
})

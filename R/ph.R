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
#' @param pi a probability vector.
#' @param T an intensity matrix.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
ph <- function(pi = c(1), T = as.matrix(- 1)) {
  if(dim(T)[1] != dim(T)[2]){stop("matrix T should be square")}
  if(length(pi) != dim(T)[1]){stop("incompatible dimensions")}
  new("ph",
      name = paste("ph(", length(pi), ")", sep = ""),
      pars = list(pi = pi, T = T),
      fit = list()
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
  t <- - rowSums(x@pars$T)
  U <- rphasetype(n, x@pars$pi, x@pars$T, t)
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
  dens <- phdensity(y, x@pars$pi, x@pars$T)
  return(cbind(y = y, dens = dens))
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
  function(x, y) {
    x@fit <- rcpp_hello_world()
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

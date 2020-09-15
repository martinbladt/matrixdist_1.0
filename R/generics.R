#' New Generic for Simulating Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x an object of the model class.
#' @param ...
#'
#' @return A realization from the matrix distribution.
#' @export
#'
#'
setGeneric("r", function(x, ...) {
  standardGeneric("r")
})

#' New Generic for the Density of Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x an object of the model class.
#' @param ...
#'
#' @return Density from the matrix distribution.
#' @export
#'
#'
setGeneric("d", function(x, ...) {
  standardGeneric("d")
})

#' New Generic for Estimating Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x an object of the model class.
#' @param y a vector of data.
#' @param ...
#'
#' @return An object of the fitted model class.
#' @export
#'
#'
setGeneric("fit", function(x, y, ...) {
  standardGeneric("fit")
})

#' New Generic for doing a likelihood ratio test between two Matrix Distribution models
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x,y objects of the model class.
#' @param ...
#'
#' @return a likelihood ratio test result.
#' @export
#'
#'
setGeneric("LRT", function(x, y, ...) {
  standardGeneric("LRT")
})

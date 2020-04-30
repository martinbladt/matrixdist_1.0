#' New Generic for Simulating Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}, \linkS4class{fph}
#'
#' @param x an object of the model class.
#' @param ...
#'
#' @return A realization from the matrix distribution.
#' @export
#'
#'
setGeneric("sim", function(x, ...) {
  standardGeneric("sim")
})

#' New Generic for Estimating Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}, \linkS4class{fph}
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
#' Methods are available for objects of class \linkS4class{ph}, \linkS4class{fph}
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

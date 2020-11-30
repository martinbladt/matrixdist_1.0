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
setGeneric("sim", function(x, ...) {
  standardGeneric("sim")
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
setGeneric("dens", function(x, ...) {
  standardGeneric("dens")
})

#' New Generic for the Hazard rate of Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x an object of the model class.
#' @param ...
#'
#' @return Hazard rate from the matrix distribution.
#' @export
#'
#'
setGeneric("haz", function(x, ...) {
  standardGeneric("haz")
})

#' New Generic for the Distribution of Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x an object of the model class.
#' @param ...
#'
#' @return CDF from the matrix distribution.
#' @export
#'
#'
setGeneric("cdf", function(x, ...) {
  standardGeneric("cdf")
})

#' New Generic for the Quantile of Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x an object of the model class.
#' @param ...
#'
#' @return Quantile from the matrix distribution.
#' @export
#'
#'
setGeneric("quan", function(x, ...) {
  standardGeneric("quan")
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

#' New Generic for Plotting Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x an object of the model class.
#' @param y a vector of data.
#' @param ...
#'
#' @export
#'
#'
setGeneric("plot", function(x, ...) {
  standardGeneric("plot")
})

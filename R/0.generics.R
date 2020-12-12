#' New Generic for Simulating Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x an object of the model class.
#' @param ... further parameters to be passed on
#'
#' @return A realization from the matrix distribution.
#'
#'
setGeneric("sim", function(x, ...) {
  standardGeneric("sim")
})

#' New Generic for Moment of Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x an object of the model class.
#' @param ... further parameters to be passed on
#'
#' @return A realization from the matrix distribution.
#'
#'
setGeneric("moment", function(x, ...) {
  standardGeneric("moment")
})


#' New Generic for Minimum of two Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x1 an object of the model class.
#' @param x2 an object of the model class.
#' @param ... further parameters to be passed on
#'
#' @return A realization from the matrix distribution.
#'
#'
setGeneric("minimum", function(x1, x2, ...) {
  standardGeneric("minimum")
})

#' New Generic for Maximum of two Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x1 an object of the model class.
#' @param x2 an object of the model class.
#' @param ... further parameters to be passed on
#'
#' @return A realization from the matrix distribution.
#'
#'
setGeneric("maximum", function(x1, x2, ...) {
  standardGeneric("maximum")
})

#' New Generic for the Density of Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x an object of the model class.
#' @param ... further parameters to be passed on
#'
#' @return Density from the matrix distribution.
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
#' @param ... further parameters to be passed on
#'
#' @return Hazard rate from the matrix distribution.
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
#' @param ... further parameters to be passed on
#'
#' @return CDF from the matrix distribution.
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
#' @param ... further parameters to be passed on
#'
#' @return Quantile from the matrix distribution.
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
#' @param ... further parameters to be passed on
#'
#' @return An object of the fitted model class.
#'
#'
setGeneric("fit", function(x, y, ...) {
  standardGeneric("fit")
})

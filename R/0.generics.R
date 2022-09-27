#' New Generic for Simulating Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return A realization from the matrix distribution.
#'
setGeneric("sim", function(x, ...) standardGeneric("sim"))

#' New Generic for Moment of Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Moment of the matrix distribution.
#'
setGeneric("moment", function(x, ...) standardGeneric("moment"))

#' New Generic for Laplace transform of Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Laplace transform of the matrix distribution.
#'
setGeneric("laplace", function(x, ...) standardGeneric("laplace"))

#' New Generic for mgf of Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Mgf transform of the matrix distribution.
#'
setGeneric("mgf", function(x, ...) standardGeneric("mgf"))

#' New Generic for pgf of Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{dph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Pgf transform of the matrix distribution.
#'
setGeneric("pgf", function(x, ...) standardGeneric("pgf"))

#' New Generic for Minimum of two Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x1 An object of the model class.
#' @param x2 An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the model class.
#'
setGeneric("minimum", function(x1, x2, ...) standardGeneric("minimum"))

#' New Generic for Maximum of two Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x1 An object of the model class.
#' @param x2 An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the model class.
#'
setGeneric("maximum", function(x1, x2, ...) standardGeneric("maximum"))

#' New Generic for Mixture of two Matrix Distributions
#'
#' Methods are available for objects of classes \linkS4class{ph} and \linkS4class{dph}.
#'
#' @param x1 An object of the model class.
#' @param x2 An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the model class.
#'
setGeneric("mixture", function(x1, x2, ...) standardGeneric("mixture"))

#' New Generic for the Density of Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Density from the matrix distribution.
#'
setGeneric("dens", function(x, ...) standardGeneric("dens"))

#' New Generic for the Hazard rate of Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Hazard rate from the matrix distribution.
#'
setGeneric("haz", function(x, ...) standardGeneric("haz"))

#' New Generic for the Distribution of Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return CDF from the matrix distribution.
#'
setGeneric("cdf", function(x, ...) standardGeneric("cdf"))

#' New Generic for the Quantile of Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Quantile from the matrix distribution.
#'
setGeneric("quan", function(x, ...) standardGeneric("quan"))

#' New Generic for Estimating Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param y A vector of data.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the fitted model class.
#'
setGeneric("fit", function(x, y, ...) standardGeneric("fit"))

#' New Generic for Regression with Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x An object of the model class.
#' @param y A vector of data.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the fitted model class.
#' @export
#'
setGeneric("reg", function(x, y, ...) standardGeneric("reg"))

#' New Generic for Evaluating Survival Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{sph}.
#'
#' @param x An object of the model class.
#' @param subject A vector of data.
#' @param ... Further parameters to be passed on.
#'
#' @export
#'
setGeneric("evaluate", function(x, subject, ...) standardGeneric("evaluate"))

#' New Generic for doing a likelihood ratio test between two Matrix Distribution
#'  models
#'
#' Methods are available for objects of class \linkS4class{ph}.
#'
#' @param x,y Objects of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return A likelihood ratio test result.
#' @export
#'
setGeneric("LRT", function(x, y, ...) standardGeneric("LRT"))

#' New Generic for obtaining the Fisher Information of Survival Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{sph}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @export
#'
setGeneric("Fisher", function(x, ...) standardGeneric("Fisher"))

#' New Generic for the transformation via rewards of a phase-type distribution
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the model class.
#'
setGeneric("TVR", function(x, ...) standardGeneric("TVR"))

#' New Generic for Regression with Matrix Distributions
#'
#' Methods are available for objects of class \linkS4class{ph}
#'
#' @param x An object of the model class.
#' @param y A vector of data.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the fitted model class.
#' @export
#'
setGeneric("MoE", function(x, y, ...) standardGeneric("MoE"))

#' New generic for the marginals of multivariate matrix distributions
#'
#' Methods are available for objects of multivariate classes.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Marginal of the matrix distribution.
#'
setGeneric("marginal", function(x, ...) standardGeneric("marginal"))

#' New generic for linear combinations of multivariate matrix distributions
#'
#' Methods are available for objects of multivariate classes.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Marginal of the matrix distribution.
#'
setGeneric("linCom", function(x, ...) standardGeneric("linCom"))

setGeneric("mean")

setGeneric("var")
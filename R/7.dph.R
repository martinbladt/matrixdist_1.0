#' Dicrete Phase Type distributions
#'
#' Class of objects for phase-type distributions.
#'
#' @slot name Name of the discrete phase-type distribution.
#' @slot pars A list comprising of the parameters.
#' @slot fit A list containing estimation information.
#'
#' @return Class object.
#' @export
#'
setClass("dph",
  slots = list(
    name = "character",
    pars = "list",
    fit = "list"
  )
)

#' Constructor Function for discrete phase-type distributions
#'
#' @param alpha A probability vector.
#' @param S A sub-transition matrix.
#' @param structure A valid dph structure ("general", "coxian", "hyperexponential", "gcoxian", "gerlang").
#' @param dimension The dimension of the dph structure (if structure is provided).
#'
#' @return An object of class \linkS4class{dph}.
#' @export
#'
#' @examples
#' dph(structure = "general", dim = 5)
#' dph(alpha = c(0.5, 0.5), S = matrix(c(0.1, 0.5, 0.5, 0.2), 2, 2))
dph <- function(alpha = NULL, S = NULL, structure = NULL, dimension = 3) {
  if (any(is.null(alpha)) & any(is.null(S)) & is.null(structure)) {
    stop("input a vector and matrix, or a structure")
  }
  if (!is.null(structure)) {
    rs <- random_structure(dimension, structure = structure)
    alpha <- rs[[1]]
    Sa <- rs[[2]]
    a <- max_diagonal(Sa * (-1)) * (1 + stats::runif(1))
    S <- (a * diag(dimension) + Sa) / a
    name <- structure
  } else {
    if (dim(S)[1] != dim(S)[2]) {
      stop("matrix S should be square")
    }
    if (length(alpha) != dim(S)[1]) {
      stop("incompatible dimensions")
    }
    name <- "custom"
  }
  methods::new("dph",
    name = paste(name, " dph(", length(alpha), ")", sep = ""),
    pars = list(alpha = alpha, S = S)
  )
}

#' Sum Method for discrete phase-type distributions
#'
#' @param e1 An object of class \linkS4class{dph}.
#' @param e2 An object of class \linkS4class{dph}.
#'
#' @return An object of class \linkS4class{dph}.
#' @export
#'
#' @examples
#' dph1 <- dph(structure = "general", dimension = 3)
#' dph2 <- dph(structure = "general", dimension = 5)
#' dph_sum <- dph1 + dph2
#' dph_sum
setMethod(
  "+", signature(e1 = "dph", e2 = "dph"),
  function(e1, e2) {
    L <- sum_dph(e1@pars$alpha, e1@pars$S, e2@pars$alpha, e2@pars$S)
    return(dph(alpha = L$alpha, S = L$S))
  }
)

#' Minimum Method for discrete phase-type distributions
#'
#' @param x1 An object of class \linkS4class{dph}.
#' @param x2 An object of class \linkS4class{dph}.
#'
#' @return An object of class \linkS4class{dph}.
#' @export
#'
#' @examples
#' dph1 <- dph(structure = "general", dimension = 3)
#' dph2 <- dph(structure = "general", dimension = 5)
#' dph_min <- minimum(dph1, dph2)
#' dph_min
setMethod(
  "minimum", signature(x1 = "dph", x2 = "dph"),
  function(x1, x2) {
    alpha <- kronecker(x1@pars$alpha, x2@pars$alpha)
    S <- kronecker(x1@pars$S, x2@pars$S)
    return(dph(alpha = alpha, S = S))
  }
)

#' Maximum Method for discrete phase-type distributions
#'
#' @param x1 An object of class \linkS4class{dph}.
#' @param x2 An object of class \linkS4class{dph}.
#'
#' @return An object of class \linkS4class{dph}.
#' @export
#'
#' @examples
#' dph1 <- dph(structure = "general", dimension = 3)
#' dph2 <- dph(structure = "general", dimension = 5)
#' dph_max <- maximum(dph1, dph2)
#' dph_max
setMethod(
  "maximum", signature(x1 = "dph", x2 = "dph"),
  function(x1, x2) {
    n1 <- length(x1@pars$alpha)
    n2 <- length(x2@pars$alpha)
    alpha <- c(kronecker(x1@pars$alpha, x2@pars$alpha), rep(0, n1 + n2))
    S1 <- rbind(kronecker(x1@pars$S, x2@pars$S), matrix(0, n1 + n2, n1 * n2))
    S2 <- rbind(kronecker(x1@pars$S, 1 - rowSums(x2@pars$S)), x1@pars$S, matrix(0, n2, n1))
    S3 <- rbind(kronecker(1 - rowSums(x1@pars$S), x2@pars$S), matrix(0, n1, n2), x2@pars$S)
    return(dph(alpha = alpha, S = cbind(S1, S2, S3)))
  }
)

#' Show Method for discrete phase-type distributions
#'
#' @param object An object of class \linkS4class{dph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "dph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@pars)
})

#' Coef Method for dph Class
#'
#' @param object An object of class \linkS4class{dph}.
#'
#' @return Parameters of dph model.
#' @export
#'
#' @examples
#' obj <- dph(structure = "general", dim = 3)
#' coef(obj)
setMethod("coef", c(object = "dph"), function(object) {
  object@pars
})


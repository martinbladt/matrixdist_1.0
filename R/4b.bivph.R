#' Bivariate phase-type distributions
#'
#' Class of objects for bivariate phase-type distributions.
#'
#' @slot name Name of the phase-type distribution.
#' @slot pars A list comprising of the parameters.
#' @slot fit A list containing estimation information.
#'
#' @return Class object.
#' @export
#'
setClass("bivph",
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

#' Constructor function for bivariate phase-type distributions
#'
#' @param alpha A probability vector.
#' @param S11 A sub-intensity matrix.
#' @param S12 A matrix.
#' @param S22 A sub-intensity matrix.
#' @param dimensions The dimensions of the bivariate phase-type (if no parameters are provided).
#'
#' @return An object of class \linkS4class{bivph}.
#' @export
#'
#' @examples
#' bivph(dimensions = c(3, 3))
#' S11 <- matrix(c(-1, .5, .5, -1), 2, 2)
#' S12 <- matrix(c(.2, .4, .3, .1), 2, 2)
#' S22 <- matrix(c(-2, 0, 1, -1), 2, 2)
#' bivph(alpha = c(.5, .5), S11, S12, S22)
bivph <- function(alpha = NULL, S11 = NULL, S12 = NULL, S22 = NULL, dimensions = c(3, 3)) {
  if (any(is.null(alpha)) & any(is.null(S11))) {
    rs <- random_structure_bivph(dimensions[1], dimensions[2])
    alpha <- rs[[1]]
    S11 <- rs[[2]]
    S12 <- rs[[3]]
    S22 <- rs[[4]]
    name <- "random"
  } else {
    if (dim(S11)[1] != dim(S11)[2]) {
      stop("matrix S11 should be square")
    } else if (dim(S22)[1] != dim(S22)[2]) {
      stop("matrix S22 should be square")
    } else if (dim(S11)[1] != dim(S12)[1]) {
      stop("incompatible dimensions of S11 and S12")
    } else if (dim(S22)[1] != dim(S12)[2]) {
      stop("incompatible dimensions of S12 and S22")
    } else if (length(alpha) != dim(S11)[1]) {
      stop("incompatible dimensions of alpha and S11")
    }
    name <- "custom"
  }
  methods::new("bivph",
    name = paste(name, " bivph(", length(alpha), ")", sep = ""),
    pars = list(alpha = alpha, S11 = S11, S12 = S12, S22 = S22)
  )
}

#' Show method for bivariate phase-type distributions
#'
#' @param object An object of class \linkS4class{bivph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "bivph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@pars)
})

#' Density method for bivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{bivph}.
#' @param y A matrix of locations.
#'
#' @return A vector containing the joint density evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' dens(obj, matrix(c(0.5, 1), ncol = 2))
setMethod("dens", c(x = "bivph"), function(x, y) {
  dens <- bivph_density(y, x@pars$alpha, x@pars$S11, x@pars$S12, x@pars$S22)
  return(dens)
})

#' Simulation method for bivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{bivph}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed bivariate
#'  phase-type vector.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' sim(obj, n = 100)
setMethod("sim", c(x = "bivph"), function(x, n = 1000) {
  p1_aux <- dim(x@pars$S11)[1]
  p2_aux <- dim(x@pars$S22)[1]
  alpha_aux <- c(x@pars$alpha, rep(0, p2_aux))
  S_aux <- merge_matrices(x@pars$S11, x@pars$S12, x@pars$S22)
  R_aux <- matrix(c(c(rep(1, p1_aux), rep(0, p2_aux)), c(rep(0, p1_aux), rep(1, p2_aux))), ncol = 2)
  U <- rMPHstar(n, alpha_aux, S_aux, R_aux)
  return(U)
})

#' Coef method for bivph class
#'
#' @param object An object of class \linkS4class{bivph}.
#'
#' @return Parameters of bivariate phase-type model.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' coef(obj)
setMethod("coef", c(object = "bivph"), function(object) {
  object@pars
})

#' Moment method for bivph class
#'
#' @param x An object of class \linkS4class{bivph}.
#' @param k A vector with the location.
#' @return An real value.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' moment(obj, c(1, 1))
setMethod("moment", c(x = "bivph"), function(x, k = c(1, 1)) {
  ee <- rep(1, nrow(x@pars$S22))
  return(factorial(k[1]) * factorial(k[2]) * x@pars$alpha %*% matrix_power(k[1] + 1, base::solve(-x@pars$S11)) %*% x@pars$S12 %*% matrix_power(k[2], base::solve(-x@pars$S22)) %*% ee)
})

#' Fit Method for bivph Class
#'
#' @param x An object of class \linkS4class{bivph}.
#' @param y A matrix with the data.
#' @param weight Vector of weights.
#' @param stepsEM Number of EM steps to be performed.
#' @param every Number of iterations between likelihood display updates.
#'
#' @return An object of class \linkS4class{bivph}.
#'
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' data <- sim(obj, n = 100)
#' fit(obj, data, stepsEM = 100, every = 50)
setMethod(
  "fit", c(x = "bivph"),
  function(x,
           y,
           weight = numeric(0),
           stepsEM = 1000,
           every = 10) {
    if (!all(y > 0)) {
      stop("data should be positive")
    }
    if (!all(weight >= 0)) {
      stop("weights should be non-negative")
    }
    if (length(weight) == 0) {
      weight <- rep(1, length(y[, 1]))
    }

    bivph_par <- x@pars
    alpha_fit <- clone_vector(bivph_par$alpha)
    S11_fit <- clone_matrix(bivph_par$S11)
    S12_fit <- clone_matrix(bivph_par$S12)
    S22_fit <- clone_matrix(bivph_par$S22)

    options(digits.secs = 4)
    cat(format(Sys.time(), format = "%H:%M:%OS"), ": EM started", sep = "")
    cat("\n", sep = "")


    for (k in 1:stepsEM) {
      EMstep_bivph(alpha_fit, S11_fit, S12_fit, S22_fit, y, weight)
      if (k %% every == 0) {
        cat("\r", "iteration:", k,
          ", logLik:", logLikelihoodbivPH(alpha_fit, S11_fit, S12_fit, S22_fit, y, weight),
          sep = " "
        )
      }
    }

    x@pars$alpha <- alpha_fit
    x@pars$S11 <- S11_fit
    x@pars$S12 <- S12_fit
    x@pars$S22 <- S22_fit

    cat("\n", format(Sys.time(), format = "%H:%M:%OS"), ": EM finalized", sep = "")
    cat("\n", sep = "")

    return(x)
  }
)

#' Marginal method for bivph class
#'
#' @param x An object of class \linkS4class{bivph}.
#' @param mar Indicator of which marginal.
#' @return An object of the of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' marginal(obj, 1)
setMethod("marginal", c(x = "bivph"), function(x, mar = 1) {
  if (mar == 1) {
    x0 <- ph(alpha = x@pars$alpha, S = x@pars$S11)
  } else {
    alpha0 <- x@pars$alpha %*% base::solve(-x@pars$S11) %*% x@pars$S12
    x0 <- ph(alpha = alpha0, S = x@pars$S22)
  }
  return(x0)
})

#' Linear Combination method for bivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{bivph}.
#' @param w A vector with non-negative entries.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' obj <- bivph(dimensions = c(3, 3))
#' linCom(obj, c(1, 0))
setMethod("linCom", c(x = "bivph"), function(x, w = c(1, 1)) {
  if (length(w) != 2) {
    stop("vector of wrong dimension")
  }
  if (any(w < 0)) {
    stop("vector with negative entries")
  }
  if (all(w == 0)) {
    stop("vector with all entries zero")
  }
  p1_aux <- dim(x@pars$S11)[1]
  p2_aux <- dim(x@pars$S22)[1]
  alpha_aux <- c(x@pars$alpha, rep(0, p2_aux))
  S_aux <- merge_matrices(x@pars$S11, x@pars$S12, x@pars$S22)
  R_aux <- matrix(c(c(rep(1, p1_aux), rep(0, p2_aux)), c(rep(0, p1_aux), rep(1, p2_aux))), ncol = 2)
  L <- linear_combination(w, alpha_aux, S_aux, R_aux)
  x0 <- ph(alpha = L$alpha, S = L$S)
  return(x0)
})
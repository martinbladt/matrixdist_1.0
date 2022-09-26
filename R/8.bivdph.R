#' Bivariate discrete phase-type distributions
#'
#' Class of objects for bivariate discrete phase-type distributions.
#'
#' @slot name Name of the discrete phase-type distribution.
#' @slot pars A list comprising of the parameters.
#' @slot fit A list containing estimation information.
#'
#' @return Class object.
#' @export
#'
setClass("bivdph",
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

#' Constructor function for bivariate discrete phase-type distributions
#'
#' @param alpha A probability vector.
#' @param S11 A sub-transition matrix.
#' @param S12 A matrix.
#' @param S22 A sub-transition matrix.
#' @param dimensions The dimensions of the bivariate discrete phase-type 
#'  (if no parameters are provided).
#'
#' @return An object of class \linkS4class{bivdph}.
#' @export
#'
#' @examples
#' bivdph(dimensions = c(3, 3))
#' S11 <- matrix(c(0.1, .5, .5, 0.1), 2, 2)
#' S12 <- matrix(c(.2, .3, .2, .1), 2, 2)
#' S22 <- matrix(c(0.2, 0, 0.1, 0.1), 2, 2)
#' bivdph(alpha = c(.5, .5), S11, S12, S22)
bivdph <- function(alpha = NULL, S11 = NULL, S12 = NULL, S22 = NULL, dimensions = c(3, 3)) {
  if (any(is.null(alpha)) & any(is.null(S11))) {
    alpha <- stats::runif(dimensions[1])
    alpha <- alpha / sum(alpha)
    aux_m <- matrix(stats::runif(dimensions[1] * (dimensions[1] + dimensions[2])), nrow = dimensions[1])
    aux_m <- aux_m / rowSums(aux_m)
    S11 <- aux_m[, 1:dimensions[1]]
    S12 <- aux_m[, -(1:dimensions[1])]
    aux_m <- matrix(stats::runif(dimensions[2] * (dimensions[2] + 1)), nrow = dimensions[2])
    aux_m <- aux_m / rowSums(aux_m)
    S22 <- aux_m[, 1:dimensions[2]]
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
  methods::new("bivdph",
    name = paste(name, " bivdph(", length(alpha), ")", sep = ""),
    pars = list(alpha = alpha, S11 = S11, S12 = S12, S22 = S22)
  )
}

#' Show method for bivariate discrete phase-type distributions
#'
#' @param object An object of class \linkS4class{bivdph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "bivdph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@pars)
})

#' Density method for bivariate discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{bivdph}.
#' @param y A matrix of locations.
#'
#' @return A vector containing the joint density evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- bivdph(dimensions = c(3, 3))
#' dens(obj, matrix(c(1, 2), ncol = 2))
setMethod("dens", c(x = "bivdph"), function(x, y) {
  dens <- bivdph_density(y, x@pars$alpha, x@pars$S11, x@pars$S12, x@pars$S22)
  return(dens)
})

#' Coef method for bivdph class
#'
#' @param object An object of class \linkS4class{bivdph}.
#'
#' @return Parameters of bivariate discrete phase-type model.
#' @export
#'
#' @examples
#' obj <- bivdph(dimensions = c(3, 3))
#' coef(obj)
setMethod("coef", c(object = "bivdph"), function(object) {
  object@pars
})

#' Simulation method for bivariate discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{bivdph}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed bivariate
#'  discrete phase-type vector.
#' @export
#'
#' @examples
#' obj <- bivdph(dimensions = c(3, 3))
#' sim(obj, n = 100)
setMethod("sim", c(x = "bivdph"), function(x, n = 1000) {
  p1_aux <- dim(x@pars$S11)[1]
  p2_aux <- dim(x@pars$S22)[1]
  alpha_aux <- c(x@pars$alpha, rep(0, p2_aux))
  S_aux <- merge_matrices(x@pars$S11, x@pars$S12, x@pars$S22)
  R_aux <- matrix(c(c(rep(1, p1_aux), rep(0, p2_aux)), c(rep(0, p1_aux), rep(1, p2_aux))), ncol = 2)
  exit_vec <- 1 - rowSums(S_aux)
  trans_mat <- cbind(S_aux, exit_vec)
  aux_vec <- rep(0, p1_aux + p2_aux + 1)
  aux_vec[p1_aux + p2_aux + 1] <- 1
  trans_mat <- rbind(trans_mat, aux_vec)
  U <- rMDPHstar(n, alpha_aux, trans_mat, R_aux)
  return(U)
})

#' Marginal method for bivdph class
#'
#' @param x An object of class \linkS4class{bivdph}.
#' @param mar Indicator of which marginal.
#' @return An object of the of class \linkS4class{dph}.
#' @export
#'
#' @examples
#' obj <- bivdph(dimensions = c(3, 3))
#' marginal(obj, 1)
setMethod("marginal", c(x = "bivdph"), function(x, mar = 1) {
  if (mar == 1) {
    x0 <- dph(alpha = x@pars$alpha, S = x@pars$S11)
  } else {
    alpha0 <- x@pars$alpha %*% base::solve(diag(ncol(x@pars$S12)) - x@pars$S11) %*% x@pars$S12
    x0 <- dph(alpha = alpha0, S = x@pars$S22)
  }
  return(x0)
})

#' Fit Method for bivdph Class
#'
#' @param x An object of class \linkS4class{bivdph}.
#' @param y A matrix with the data.
#' @param weight Vector of weights.
#' @param stepsEM Number of EM steps to be performed.
#' @param every Number of iterations between likelihood display updates.
#'
#' @return An object of class \linkS4class{bivdph}.
#'
#' @export
#'
#' @examples
#' obj <- bivdph(dimensions = c(3, 3))
#' data <- sim(obj, n = 100)
#' fit(obj, data, stepsEM = 100, every = 50)
setMethod(
  "fit", c(x = "bivdph"),
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

    bivdph_par <- x@pars
    alpha_fit <- clone_vector(bivdph_par$alpha)
    S11_fit <- clone_matrix(bivdph_par$S11)
    S12_fit <- clone_matrix(bivdph_par$S12)
    S22_fit <- clone_matrix(bivdph_par$S22)

    options(digits.secs = 4)
    cat(format(Sys.time(), format = "%H:%M:%OS"), ": EM started", sep = "")
    cat("\n", sep = "")

    for (k in 1:stepsEM) {
      EMstep_bivdph(alpha_fit, S11_fit, S12_fit, S22_fit, y, weight)
      if (k %% every == 0) {
        cat("\r", "iteration:", k,
          ", logLik:", logLikelihoodbivDPH(alpha_fit, S11_fit, S12_fit, S22_fit, y, weight),
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

#' MoE Method for bivdph Class
#'
#' @param x An object of class \linkS4class{bivdph}.
#' @param formula A regression formula.
#' @param y A matrix of observations.
#' @param data A data frame of covariates.
#' @param alpha_vecs Matrix of initial probabilities.
#' @param weight Vector of weights.
#' @param stepsEM Number of EM steps to be performed.
#' @param every Number of iterations between likelihood display updates.
#' @param rand_init Random initiation in the R-step.
#'
#' @return An object of class \linkS4class{sph}.
#'
#' @export
#'
setMethod(
  "MoE", c(x = "bivdph"),
  function(x,
           formula,
           y,
           data,
           alpha_vecs = NULL,
           weight = numeric(0),
           stepsEM = 1000,
           every = 10,
           rand_init = TRUE) {
    p <- length(x@pars$alpha)
    frame <- stats::model.frame(formula, data = data)
    n <- nrow(frame)
    d <- ncol(frame) - 1
    if (is.null(alpha_vecs)) alpha_vecs <- matrix(x@pars$alpha, ncol = p, nrow = n, byrow = TRUE)
    if (length(weight) == 0) weight <- rep(1, n)
    S11_fit <- clone_matrix(x@pars$S11)
    S12_fit <- clone_matrix(x@pars$S12)
    S22_fit <- clone_matrix(x@pars$S22)
    c <- c()
    for (i in 1:p) c <- c(c, rep(i, n))
    extended_x <- matrix(t(as.matrix(frame[, -1])), nrow = n * p, ncol = d, byrow = TRUE)
    dm <- data.frame(Class = c, extended_x)
    names(dm)[-1] <- names(frame)[-1]
    ndm <- data.frame(dm[dm$Class == 1, -1])
    names(ndm) <- names(dm)[-1]
    for (k in 1:stepsEM) {
      B_matrix <- EMstep_bivdph_MoE(alpha_vecs, S11_fit, S12_fit, S22_fit, y, weight)
      wt <- reshape2::melt(B_matrix)[, 3]
      wt[wt < 1e-22] <- wt[wt < 1e-22] + 1e-22
      if (k == 1 | rand_init == TRUE) {
        multinom_model <- nnet::multinom(Class ~ ., data = dm, weights = wt, trace = F)
      } else {
        multinom_model <- nnet::multinom(Class ~ ., data = dm, weights = wt, trace = F, Wts = multinom_model$wts)
      }
      alpha_vecs <- stats::predict(multinom_model, type = "probs", newdata = ndm)
      if (k %% every == 0) {
        ll <- logLikelihoodbivDPH_MoE(alpha_vecs, S11_fit, S12_fit, S22_fit, y, weight)
        cat("\r", "iteration:", k, ", logLik:", ll, sep = " ")
      }
    }
    cat("\n", sep = "")
    return(list(alpha = alpha_vecs, S11 = S11_fit, S12 = S11_fit, S12 = S11_fit, mm = multinom_model))
  }
)


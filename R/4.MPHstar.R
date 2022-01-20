
#' Multivariate Phase Type distributions obtained by transformation via rewards
#'
#' Class of objects for multivariate phase type distributions.
#'
#' @slot name Name of the phase type distribution.
#' @slot pars A list comprising of the parameters.
#'
#' @export
#'
setClass("MPHstar",
  slots = list(
    name = "character",
    pars = "list"
  )
)

#' Constructor Function for multivariate phase type distributions (MPH* class)
#'
#' @param alpha A probability vector.
#' @param S A sub-intensity matrix.
#' @param R A compatible (non-negative) reward matrix.
#' @param structure A valid ph structure.
#' @param dimension The dimension of the ph structure (if provided).
#' @param variables The number of desired marginals.
#'
#' @return An object of class \linkS4class{MPHstar}.
#' @export
#'
MPHstar <- function(alpha = NULL,
                    S = NULL,
                    structure = NULL,
                    dimension = 3,
                    R = NULL,
                    variables = 2) {
  if (any(is.null(alpha)) & any(is.null(S)) & is.null(structure)) {
    stop("input a vector and matrix, or a structure")
  }

  if (!all(is.null(R))) {
    if (R %*% rep(1, ncol(R)) != rep(1, ncol(R))) stop("Rows of the reward matrix must sum to 1")
    Rname <- "custom"
  }

  if (!is.null(structure)) {
    rs <- random_structure(p = dimension, structure = structure, scale_factor = 1)
    alpha <- rs[[1]]
    S <- rs[[2]]
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

  if (all(is.null(R))) {
    R <- random_reward(dim(S)[1], variables)
    Rname <- "random"
  } else {
    if (dim(R)[1] != dim(S)[1]) {
      stop("The reward matrix R has wrong dimensions: number of rows does not match ph dimension")
      Rname <- "custom"
    }
  }


  methods::new("MPHstar",
    name = paste("mph based on a ", name, " ph and a ", Rname, " reward matrix", sep = ""),
    pars = list(alpha = alpha, S = S, R = R)
  )
}

#' Show Method for multivariate phase type distributions
#'
#' @param object An object of class \linkS4class{MPHstar}.
#' @export
#'
setMethod("show", "MPHstar", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("number of marginal distributions: ", dim(object@pars$R)[2], "\n", sep = "")
})

#' Simulation Method for multivariate phase type distributions
#'
#' @param x An object of class \linkS4class{MPHstar}.
#' @param n Desired sample size for each marginal.
#'
#' @return A matrix of sample data for each marginal.
#' @export
#'
setMethod("sim", c(x = "MPHstar"), function(x, n = 1000) {
  U <- rMPHstar(n, x@pars$alpha, x@pars$S, x@pars$R)
  return(U)
})

#' Transformation Via Rewards method for multivariate phase type distributions
#'
#' @param x An object of class \linkS4class{MPHstar}.
#'
#' @return A list with marginal distributions, obtained by transformation via rewards.
#' @export
#'
setMethod("TVR", c(x = "MPHstar"), function(x) {
  alpha <- x@pars$alpha
  S <- x@pars$S
  R <- x@pars$R
  Q <- embedded_mc(S)

  marginal <- transf_via_rew(R, Q, alpha, S)
  return(marginal)
})

#' Find weight of observations
#'
#' @param x A vector of observations from which we want to know their weights.
#'
#' @return A matrix with unique observations as fisrt column and associated weights for second column.
#' @export
#'
find_weight <- function(x) {
  cum_weight <- numeric(0)
  x <- sort(x)
  unique <- unique(x)
  for (i in unique) {
    if (i > 0) {
      cum_weight <- c(cum_weight, length(which(x == i)))
    }
  }

  return(cbind(unique, cum_weight, deparse.level = 0))
}

#' Prepare data for the MPHstar_EMstep_UNI
#'
#' @param y A matrix with marginal observations, each column corresponds to a marginal.
#' @param w A matrix of weights, each column corresponds to a marginal.
#' @param rc A matrix with indication if an observation is right-censored (0 if right-censored ).
#'
#' @return For summed and marginal observations we have a list with matrices of unique observations and their associated weights, separated by uncensored and right-censored data.
#' @export
#'
MPHstar_data_aggregation <- function(y, w = numeric(0), rc = numeric(0)) {
  mat <- list()

  if (is.matrix(y) & length(y) > 1) {
    if (length(rc) == 0) {
      rc <- matrix(rep(1, nrow(y) * ncol(y)), nrow(y), ncol(y))
    } # no right censoring

    sumData <- rowSums(y)
    un_obs <- which(rowSums(rc) == ncol(y)) # all summed observations that are uncensored
    rc_obs <- which(rowSums(rc) < ncol(y)) # all summed observations that right-censored

    n1 <- length(un_obs)
    n2 <- length(rc_obs)

    if (n1 > 1 & n2 > 1) {
      mat[[1]] <- list(un = find_weight(sumData[un_obs]), rc = find_weight(sumData[rc_obs]))
    } else if (n1 == 1 & n2 > 1) {
      mat[[1]] <- list(un = find_weight(sum(y[, un_obs])), rc = find_weight(sumData[rc_obs]))
    } else if (n1 == 0 & n2 > 1) {
      mat[[1]] <- list(un = find_weight(NULL), rc = find_weight(sumData[rc_obs]))
    } else if (n1 > 1 & n2 == 1) {
      mat[[1]] <- list(un = find_weight(sumData[un_obs]), rc = find_weight(sum(y[, rc_obs])))
    } else if (n1 > 1 & n2 == 0) {
      mat[[1]] <- list(un = find_weight(sumData[un_obs]), rc = find_weight(NULL))
    } else if (n1 == 1 & n2 == 1) {
      mat[[1]] <- list(un = find_weight(sum(y[, un_obs])), rc = find_weight(sum(y[, rc_obs])))
    }
  } else if (is.matrix(y) & length(y) == 1) {
    stop("Please input a matrix of observations")
  }

  for (i in 1:ncol(y)) {
    m <- 1 + i
    m_y <- y[, i]
    m_rc <- rc[, i]

    if (length(w) == 0) {
      mat[[m]] <- list(
        un = find_weight(m_y[which(m_rc == 1)]),
        rc = find_weight(m_y[which(m_rc == 0)])
      )
    }
    if (length(w) > 0) {
      m_w <- w[, i]

      mat[[m]] <- list(
        un = cbind(m_y[which(m_rc == 1)], m_w[which(m_rc == 1)]),
        rc = cbind(m_y[which(m_rc == 0)], m_w[which(m_rc == 0)])
      )
    }
  }

  return(mat)
}

#' Fit Method for mph Class
#'
#' @param x An object of class \linkS4class{MPHstar}.
#' @param y A matrix of marginal data.
#' @param rc A matrix with indication about right-censoring of an observation. 0 for uncensored and 1 for right-censored.
#' @param weight A matrix of marginal weights.
#' @param stepsEM The number of EM steps to be performed, defaults to 1000.
#' @param uni_epsilon The epsilon parameter for the uniformization method, defaults to 1e-4.
#' @param zero_tol The smallest value that a reward can take (to avoid numerical instability), defaults to 1e-4.
#' @param every The number of iterations between likelihood display updates. The originating distribution is used, given that there is no explicit density.
#' @param plot Boolean that determines if the plot of the loglikelihood evolution is plotted, defaults to False.
#' @param r The sub-sampling proportion for stochastic EM, defaults to 1.
#' @param replace Boolean that determines if sub-sampling is done with replacement or not, defaults to False.
#'
#' @return An object of class \linkS4class{MPHstar}.
#'
#' @importFrom grDevices dev.off
#' @importFrom graphics hist legend lines
#' @importFrom utils head
#'
#' @export
#'
setMethod(
  "fit", c(x = "MPHstar", y = "ANY"),
  function(x,
           y,
           rc = numeric(0),
           weight = numeric(0),
           stepsEM = 1000,
           uni_epsilon = 1e-4,
           zero_tol = 1e-4,
           every = 100,
           plot = F,
           r = 1,
           replace = F) {
    if (!is.matrix(y)) {
      stop("Please use a matrix with marginal observations")
    }
    if (!all(y > 0)) {
      stop("data should be positive")
    }
    if (!all(weight >= 0)) {
      stop("weights should be non-negative")
    }
    if (!all(x@pars$R >= 0)) {
      stop("rewards should be non-negative")
    }
    if (all(x@pars$R == 0)) {
      stop("A reward matrix of zeros results in an atom at zero with probability 100%")
    }
    if (r <= 0 && r > 1) {
      stop("sub-sampling proportion is invalid, please input a r in (0,1]")
    }

    is_mph <- methods::is(x, "MPHstar")

    A <- MPHstar_data_aggregation(y, weight, rc)
    B <- A
    C <- A

    ph_par <- x@pars
    alpha_fit <- clone_vector(ph_par$alpha)
    S_fit <- clone_matrix(ph_par$S)
    R_fit <- clone_matrix(ph_par$R)

    log_check <- rep(0, stepsEM / every)

    if (is_mph) {
      options(digits.secs = 4)
      cat(format(Sys.time(), format = "%H:%M:%OS"), ": EM started", sep = "")
      cat("\n", sep = "")

      for (k in 1:stepsEM) {
        if (r < 1) {
          r_rc <- numeric(0)
          r_weight <- numeric(0)

          indices <- sample(1:nrow(y), size = floor(r * nrow(y)), replace = replace)
          r_y <- y[indices, ]

          if (length(rc) > 0) {
            r_rc <- rc[indices, ]
          }
          if (length(weight) > 0) {
            r_weight <- weight[indices, ]
          }

          B <- MPHstar_data_aggregation(r_y, r_weight, r_rc)
        }

        epsilon <- if (!is.na(uni_epsilon)) {
          uni_epsilon
        } else {
          1e-4
        }

        MPHstar_EMstep_UNI(epsilon, zero_tol, alpha_fit, S_fit, R_fit, B) # performs a EM step, changes alpha, S and R

        if (length(B[[1]]$rc) > 0) {
          log_lik <- matrixdist:::logLikelihoodPH_UNI(epsilon, alpha_fit, S_fit, C[[1]]$un[1, ], C[[1]]$un[2, ], C[[1]]$rc[1, ], C[[1]]$rc[2, ])
        } else {
          log_lik <- matrixdist:::logLikelihoodPH_UNI(epsilon, alpha_fit, S_fit, C[[1]]$un[1, ], C[[1]]$un[2, ], numeric(0), numeric(0))
        }

        log_check[k] <- log_lik

        if (k %% every == 0) {
          if (length(B[[1]]$rc) > 0) {
            cat("\r", "iteration:", k,
              ", logLik:", log_lik,
              sep = " "
            )
          } else {
            cat("\r", "iteration:", k,
              ", logLik:", log_lik,
              sep = " "
            )
          }
        }
      }
    }
    if (plot == T) {
      plot(log_check, type = "l", main = "Evolution of incomplete loglikelihood", ylab = "loglik", xlab = "step EM")
    }

    cat("\n", sep = "")
    x@pars$alpha <- alpha_fit
    x@pars$S <- S_fit
    x@pars$R <- R_fit

    cat("\n", format(Sys.time(), format = "%H:%M:%OS"), ": EM finalized", sep = "")
    cat("\n", sep = "")

    return(x)
  }
)

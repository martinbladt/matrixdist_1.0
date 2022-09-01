#' MoE Method for ph Class
#'
#' @param x An object of class \linkS4class{ph}.
#' @param formula A regression formula.
#' @param data A data frame.
#' @param inhom A list with the inhomogeneity functions.
#' @param alpha_vecs <atrix of initial probabilities.s
#' @param weight Vector of weights.
#' @param delta Right-censoring indicator.
#' @param stepsEM Number of EM steps to be performed.
#' @param optim_method Method to use in gradient optimization
#' @param maxit Maximum number of iterations when optimizing g function.
#' @param reltol Relative tolerance when optimizing g function.
#' @param every Number of iterations between likelihood display updates.
#' @param rand_init Random initiation in the R-step.
#'
#' @return An object of class \linkS4class{sph}.
#'
#' @importFrom methods is new
#' @importFrom stats optim
#' @importFrom utils tail
#'
#' @export
#'
setMethod(
  "MoE", c(x = "ph"),
  function(x,
           formula,
           data,
           inhom = NULL,
           alpha_vecs = NULL,
           weight = numeric(0),
           delta = numeric(0),
           stepsEM = 1000,
           optim_method = "BFGS",
           maxit = 50,
           reltol = 1e-8,
           every = 10,
           rand_init = TRUE) {
    control <- if (optim_method == "BFGS") {
      list(
        maxit = maxit,
        factr = reltol
      )
    } else {
      list(
        maxit = maxit,
        reltol = reltol
      )
    }
    inh <- !is.null(inhom)
    if (inh) {
      g_inv <- inhom$g_inv
      lambda <- inhom$lambda
      theta <- inhom$theta
      mLL <- function(theta, g_inv, lambda, alpha1, alpha2, S, y, w, yc, wc) {
        return(-logLikelihoodPH_MoE(alpha1, alpha2, S, g_inv(y, theta), w, g_inv(yc, theta), wc) -
          sum(w * log(lambda(y, theta))))
      }
    }
    p <- length(x@pars$alpha)
    frame <- stats::model.frame(formula, data = data)
    n <- nrow(frame)
    d <- ncol(frame) - 1
    if (is.null(alpha_vecs)) alpha_vecs <- matrix(x@pars$alpha, ncol = p, nrow = n, byrow = TRUE)
    if (length(weight) == 0) weight <- rep(1, n)
    if (length(delta) == 0) delta <- rep(1, n)
    S_fit <- clone_matrix(x@pars$S)
    c <- c()
    for (i in 1:p) c <- c(c, rep(i, n)) # classes for the B matrix observations
    extended_x <- matrix(t(as.matrix(frame[, -1])), nrow = n * p, ncol = d, byrow = TRUE) # extended form of covariates
    dm <- data.frame(Class = c, extended_x)
    names(dm)[-1] <- names(frame)[-1]
    ndm <- data.frame(dm[dm$Class == 1, -1])
    names(ndm) <- names(dm)[-1]
    for (k in 1:stepsEM) {
      if (inh) {
        B_matrix <- EMstep_MoE_PADE(alpha_vecs, S_fit, g_inv(x = frame[delta == 1, 1], theta = theta), weight[delta == 1], g_inv(x = frame[delta == 0, 1], theta = theta), weight[delta == 0])[[1]]
        wt <- reshape2::melt(B_matrix)[, 3]
        wt[wt < 1e-22] <- wt[wt < 1e-22] + 1e-22
        if (k == 1 | rand_init == TRUE) {
          multinom_model <- nnet::multinom(Class ~ ., data = dm, weights = wt, trace = F)
        } else {
          multinom_model <- nnet::multinom(Class ~ ., data = dm, weights = wt, trace = F, Wts = multinom_model$wts)
        }
        alpha_vecs <- stats::predict(multinom_model, type = "probs", newdata = ndm)
        opt <- suppressWarnings(optim(
          par = theta, fn = mLL,
          g_inv = g_inv,
          lambda = lambda,
          alpha1 = alpha_vecs[delta == 1, ],
          alpha2 = alpha_vecs[delta == 0, ],
          S = S_fit,
          y = frame[delta == 1, 1],
          w = weight[delta == 1],
          yc = frame[delta == 0, 1],
          wc = weight[delta == 0],
          method = optim_method,
          control = control
        ))
        theta <- opt$par
        if (k %% every == 0) {
          cat("\r", "iteration:", k, ", logLik:", -opt$value, sep = " ")
        }
      } else {
        B_matrix <- EMstep_MoE_PADE(alpha_vecs, S_fit, frame[delta == 1, 1], weight[delta == 1], frame[delta == 0, 1], weight[delta == 0])[[1]]
        wt <- reshape2::melt(B_matrix)[, 3]
        wt[wt < 1e-22] <- wt[wt < 1e-22] + 1e-22
        if (k == 1 | rand_init == TRUE) {
          multinom_model <- nnet::multinom(Class ~ ., data = dm, weights = wt, trace = F)
        } else {
          multinom_model <- nnet::multinom(Class ~ ., data = dm, weights = wt, trace = F, Wts = multinom_model$wts)
        }
        alpha_vecs <- stats::predict(multinom_model, type = "probs", newdata = ndm)
        if (k %% every == 0) {
          ll <- logLikelihoodPH_MoE(alpha_vecs[delta == 1, ], alpha_vecs[delta == 0, ], S_fit, frame[delta == 1, 1], weight[delta == 1], frame[delta == 0, 1], weight[delta == 0])
          cat("\r", "iteration:", k, ", logLik:", ll, sep = " ")
        }
      }
    }
    if (inh) inhom$theta <- theta
    cat("\n", sep = "")
    return(list(alpha = alpha_vecs, S = S_fit, mm = multinom_model, inhom = inhom))
  }
)

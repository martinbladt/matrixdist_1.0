#' Phase-type distributions
#'
#' Class of objects for phase-type distributions.
#'
#' @slot name Name of the phase-type distribution.
#' @slot pars A list comprising of the parameters.
#' @slot fit A list containing estimation information.
#'
#' @return Class object.
#' @export
#'
setClass("ph",
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

.normalize_em_methods <- function(methods) {
  if (length(methods) == 1L) {
    methods <- rep(methods, 2L)
  }
  if (length(methods) != 2L) {
    stop("methods should have length 1 or 2")
  }
  methods <- toupper(trimws(as.character(methods)))
  allowed <- c("RK", "UNI", "PADE")
  if (any(!methods %in% allowed)) {
    stop("methods should be one of: RK, UNI, PADE")
  }
  methods
}

#' Constructor function for phase-type distributions
#'
#' @param alpha A probability vector.
#' @param S A sub-intensity matrix.
#' @param structure A valid ph structure: `"general"`, `"coxian"`,
#' `"hyperexponential"`, `"gcoxian"`, `"gerlang"`, `"erlang"`, or `"merlang"`.
#' @param dimension The dimension of the ph structure (if structure is provided).
#' @param block_sizes Optional integer vector with Erlang block sizes when
#'  `structure = "merlang"` (or aliases `"mixederlang"`, `"mixed_erlang"`).
#' @param probs Optional mixture weights used only for `structure = "merlang"`.
#' @param rates Optional Erlang rates used only for `structure = "merlang"`.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' ph(structure = "gcoxian", dimension = 5)
#' ph(structure = "merlang", block_sizes = c(2, 3))
#' ph(alpha = c(.5, .5), S = matrix(c(-1, .5, .5, -1), 2, 2))
ph <- function(alpha = NULL, S = NULL, structure = NULL, dimension = 3,
               block_sizes = NULL, probs = NULL, rates = NULL) {
  has_structure <- !is.null(structure)
  has_alpha <- !is.null(alpha)
  has_S <- !is.null(S)
  
  if (!has_structure && (!has_alpha || !has_S)) {
    stop("input a vector and matrix, or a structure")
  }
  
  if (has_structure) {
    if (!is.character(structure) || length(structure) != 1L) {
      stop("structure should be a character string of length one")
    }
    st <- tolower(structure)
    if (st %in% c("merlang", "mixederlang", "mixed_erlang")) {
      return(merlang_ph(block_sizes = block_sizes, probs = probs, rates = rates))
    }
    rs <- random_structure(dimension, structure = st)
    alpha <- rs[[1]]
    S <- rs[[2]]
    name <- st
  } else {
    if (!is.matrix(S)) {
      stop("S should be a matrix")
    }
    if (nrow(S) != ncol(S)) {
      stop("matrix S should be square")
    }
    if (length(alpha) != nrow(S)) {
      stop("incompatible dimensions")
    }
    name <- "custom"
  }
  
  methods::new(
    "ph",
    name = paste(name, " ph(", length(alpha), ")", sep = ""),
    pars = list(alpha = alpha, S = S)
  )
}

#' Sum method for phase-type distributions
#'
#' @param e1 An object of class \linkS4class{ph}.
#' @param e2 An object of class \linkS4class{ph}.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' ph1 <- ph(structure = "general", dimension = 3)
#' ph2 <- ph(structure = "gcoxian", dimension = 5)
#' ph_sum <- ph1 + ph2
#' ph_sum
setMethod(
  "+", signature(e1 = "ph", e2 = "ph"),
  function(e1, e2) {
    if (methods::is(e1, "iph") | methods::is(e2, "iph")) {
      stop("objects to be added should be ph")
    }
    L <- sum_ph(e1@pars$alpha, e1@pars$S, e2@pars$alpha, e2@pars$S)
    ph(alpha = L$alpha, S = L$S)
  }
)

kronecker_sum <- function(A, B) {
  n <- nrow(A)
  m <- nrow(B)
  kronecker(A, diag(m)) + kronecker(diag(n), B)
}

#' Minimum method for phase-type distributions
#'
#' @param x1 An object of class \linkS4class{ph}.
#' @param x2 An object of class \linkS4class{ph}.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' ph1 <- ph(structure = "general", dimension = 3)
#' ph2 <- ph(structure = "gcoxian", dimension = 5)
#' ph_min <- minimum(ph1, ph2)
#' ph_min
setMethod(
  "minimum", signature(x1 = "ph", x2 = "ph"),
  function(x1, x2) {
    alpha <- kronecker(x1@pars$alpha, x2@pars$alpha)
    S <- kronecker_sum(x1@pars$S, x2@pars$S)
    ph(alpha = alpha, S = S)
  }
)

#' Maximum method for phase-type distributions
#'
#' @param x1 An object of class \linkS4class{ph}.
#' @param x2 An object of class \linkS4class{ph}.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' ph1 <- ph(structure = "general", dimension = 3)
#' ph2 <- ph(structure = "gcoxian", dimension = 5)
#' ph_max <- maximum(ph1, ph2)
#' ph_max
setMethod(
  "maximum", signature(x1 = "ph", x2 = "ph"),
  function(x1, x2) {
    n1 <- length(x1@pars$alpha)
    n2 <- length(x2@pars$alpha)
    alpha <- c(kronecker(x1@pars$alpha, x2@pars$alpha), rep(0, n1 + n2))
    S1 <- rbind(kronecker_sum(x1@pars$S, x2@pars$S), matrix(0, n1 + n2, n1 * n2))
    S2 <- rbind(kronecker(diag(n1), -rowSums(x2@pars$S)), x1@pars$S, matrix(0, n2, n1))
    S3 <- rbind(kronecker(-rowSums(x1@pars$S), diag(n2)), matrix(0, n1, n2), x2@pars$S)
    ph(alpha = alpha, S = cbind(S1, S2, S3))
  }
)

#' Mixture method for phase-type distributions
#'
#' @param x1 An object of class \linkS4class{ph}.
#' @param x2 An object of class \linkS4class{ph}.
#' @param prob Probability for first object.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' ph1 <- ph(structure = "general", dimension = 3)
#' ph2 <- ph(structure = "gcoxian", dimension = 5)
#' ph_mix <- mixture(ph1, ph2, 0.5)
#' ph_mix
setMethod(
  "mixture", signature(x1 = "ph", x2 = "ph"),
  function(x1, x2, prob) {
    n1 <- length(x1@pars$alpha)
    n2 <- length(x2@pars$alpha)
    alpha <- c(prob * x1@pars$alpha, (1 - prob) * x2@pars$alpha)
    S1 <- rbind(x1@pars$S, matrix(0, n2, n1))
    S2 <- rbind(matrix(0, n1, n2), x2@pars$S)
    ph(alpha = alpha, S = cbind(S1, S2))
  }
)

#' Moment method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param k A positive integer (moment order).
#'
#' @return The raw moment of the \linkS4class{ph} (or underlying \linkS4class{ph}) object.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- ph(structure = "general", dimension = 3)
#' moment(obj, 2)
setMethod(
  "moment", signature(x = "ph"),
  function(x, k = 1) {
    if (k <= 0) {
      stop("k should be positive")
    }
    if ((k %% 1) != 0) {
      stop("k should be an integer")
    }
    if (methods::is(x, "iph")) {
      warning("moment of undelying ph structure is provided for iph objects")
    }
    prod <- matrix_power(k, solve(-x@pars$S))
    factorial(k) * sum(x@pars$alpha %*% prod)
  }
)

#' Mean method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#'
#' @return The raw first moment of the \linkS4class{ph} (or underlying \linkS4class{ph}) object.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- ph(structure = "general", dimension = 3)
#' mean(obj)
setMethod(
  "mean", signature(x = "ph"),
  function(x) {
    if (methods::is(x, "iph")) {
      warning("moment of undelying ph structure is provided for iph objects")
    }
    m <- solve(-x@pars$S)
    sum(x@pars$alpha %*% m)
  }
)

#' Var method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#'
#' @return The variance of the \linkS4class{ph} (or underlying \linkS4class{ph}) object.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- ph(structure = "general", dimension = 3)
#' var(obj)
setMethod(
  "var", signature(x = "ph"),
  function(x) {
    if (methods::is(x, "iph")) {
      warning("variance of undelying ph structure is provided for iph objects")
    }
    m <- solve(-x@pars$S)
    m2 <- matrix_power(2, m)
    2 * sum(x@pars$alpha %*% m2) - (sum(x@pars$alpha %*% m))^2
  }
)

#' Laplace method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param r A vector of real values.
#'
#' @return The Laplace transform of the \linkS4class{ph}
#'  (or underlying \linkS4class{ph}) object at the given locations.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- ph(structure = "general", dimension = 3)
#' laplace(obj, 3)
setMethod(
  "laplace", signature(x = "ph"),
  function(x, r) {
    if (methods::is(x, "iph")) {
      warning("Laplace transform of undelying ph structure is provided for iph objects")
    }
    lim <- max(Re(eigen(x@pars$S)$values))
    if (any(r <= lim)) {
      stop("r should be above the largest real eigenvalue of S")
    }
    ph_laplace(r, x@pars$alpha, x@pars$S)
  }
)

#' Mgf method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param r A vector of real values.
#'
#' @return The mgf of the \linkS4class{ph} (or underlying \linkS4class{ph}) object
#'  at the given locations.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- ph(structure = "general", dimension = 3)
#' mgf(obj, 0.4)
setMethod(
  "mgf", signature(x = "ph"),
  function(x, r) {
    if (methods::is(x, "iph")) {
      warning("mgf of undelying ph structure is provided for iph objects")
    }
    lim <- -max(Re(eigen(x@pars$S)$values))
    if (any(r > lim)) {
      stop("r should be below the negative largest real eigenvalue of S")
    }
    ph_laplace(-r, x@pars$alpha, x@pars$S)
  }
)

#' Show method for phase-type distributions
#'
#' @param object An object of class \linkS4class{ph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "ph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@pars)
})

#' Simulation method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed phase-type variables.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' sim(obj, n = 100)
setMethod("sim", c(x = "ph"), function(x, n = 1000) {
  rphasetype(n, x@pars$alpha, x@pars$S)
})

#' Density method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param y A vector of locations.
#'
#' @return A vector containing the density evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' dens(obj, c(1, 2, 3))
setMethod("dens", c(x = "ph"), function(x, y) {
  y_inf <- (y == Inf)
  dens <- y
  dens[!y_inf] <- phdensity(y, x@pars$alpha, x@pars$S)
  dens[y_inf] <- 0
  dens
})

#' Distribution method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param q A vector of locations.
#' @param lower.tail Logical parameter specifying whether lower tail (CDF) or
#' upper tail is computed.
#'
#' @return A vector containing the CDF evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' cdf(obj, c(1, 2, 3))
setMethod("cdf", c(x = "ph"), function(x,
                                       q,
                                       lower.tail = TRUE) {
  q_inf <- (q == Inf)
  cdf <- q
  cdf[!q_inf] <- phcdf(q[!q_inf], x@pars$alpha, x@pars$S, lower.tail)
  cdf[q_inf] <- as.numeric(1 * lower.tail)
  cdf
})

#' Hazard rate method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param y A vector of locations.
#'
#' @return A vector containing the hazard rate evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' haz(obj, c(1, 2, 3))
setMethod("haz", c(x = "ph"), function(x, y) {
  d <- dens(x, y)
  s <- cdf(x, y, lower.tail = FALSE)
  d / s
})

#' Quantile method for phase-type distributions
#'
#' @param x An object of class \linkS4class{ph}.
#' @param p A vector of probabilities.
#'
#' @return A vector containing the quantile evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' quan(obj, c(0.5, 0.9, 0.99))
setMethod("quan", c(x = "ph"), function(x,
                                        p) {
  quan <- numeric(length(p))
  for (i in seq_along(p)) {
    quan[i] <- stats::uniroot(f = function(q) p[i] - cdf(x, 1 / (1 - q) - 1), interval = c(0, 1))$root
  }
  1 / (1 - quan) - 1
})

#' Fit method for ph class
#'
#' @param x An object of class \linkS4class{ph}.
#' @param y Vector or data.
#' @param weight Vector of weights.
#' @param rcen Vector of right-censored observations. If a matrix is passed as argument, then observations are considered interval-censored
#' @param rcenweight Vector of weights for right-censored observations.
#' @param stepsEM Number of EM steps to be performed.
#' @param methods Methods to use for matrix exponential calculation: RM, UNI or PADE.
#' @param rkstep Runge-Kutta step size (optional).
#' @param uni_epsilon Epsilon parameter for uniformization method.
#' @param optim_method Optimization method for optim function.
#' @param maxit Maximum number of iterations when optimizing g function.
#' @param reltol Relative tolerance when optimizing g function.
#' @param every Number of iterations between likelihood display updates.
#' @param r Sub-sampling proportion for stochastic EM, defaults to 1.
#' @param burn Percentage of iterations (from 0 to 100) run on histogram
#'  midpoints of `y` with bin counts as weights, before switching to the
#'  standard fit on original data. Available only for uncensored data.
#' @param erlang Logical flag for exact Erlang EM updates with one repeated rate.
#'  If `NULL`, this is auto-detected from the model name.
#' @param merlang_blocks Optional integer vector with Erlang block sizes for
#'  exact mixture-of-Erlangs EM updates. If `NULL`, it is auto-read from object
#'  attributes when available.
#' @param plot_progress Logical flag. If `TRUE`, a progress plot is redrawn
#'  during EM: histogram plus fitted density for uncensored data, or Kaplan-Meier
#'  plus fitted survival for right-censored data.
#' @param plot_every Positive integer. Plot refresh frequency in iterations when
#'  `plot_progress = TRUE`. Defaults to 100.
#' @param plot_breaks Number of histogram breaks used in progress plotting.
#' @param plot_pause Small pause (in seconds) after each plot refresh to let
#'  interactive devices repaint in real time.
#'
#' @return An object of class \linkS4class{ph}.
#'
#' @importFrom grDevices dev.off
#' @importFrom graphics hist legend lines
#' @importFrom utils head
#'
#' @export
#'
#' @examples
#' obj <- iph(ph(structure = "general", dimension = 2), gfun = "weibull", gfun_pars = 2)
#' data <- sim(obj, n = 100)
#' fit(obj, data, stepsEM = 100, every = 20)
setMethod(
  "fit", c(x = "ph", y = "ANY"),
  function(x,
           y,
           weight = numeric(0),
           rcen = numeric(0),
           rcenweight = numeric(0),
           stepsEM = 1000,
           methods = c("RK", "RK"),
           rkstep = NA,
           uni_epsilon = NA,
           maxit = 100,
           reltol = 1e-8,
           every = 100,
           optim_method = "BFGS",
           r = 1,
           burn = 0,
           erlang = NULL,
           merlang_blocks = NULL,
           plot_progress = FALSE,
           plot_every = 100,
           plot_breaks = 50,
           plot_pause = 0.01) {
    methods <- .normalize_em_methods(methods)
    control <- if (optim_method == "BFGS") {
      list(
        maxit = maxit,
        factr = reltol,
        fnscale = -1
      )
    } else {
      list(
        maxit = maxit,
        reltol = reltol,
        fnscale = -1
      )
    }
    rightCensored <- is.vector(rcen)
    has_censoring <- if (is.matrix(rcen)) nrow(rcen) > 0 else length(rcen) > 0
    if(rightCensored){
      EMstep <- eval(parse(text = paste("EMstep_", methods[1], sep = "")))
    }else if(is.matrix(rcen)){
      EMstep <- eval(parse(text = "EMstep_UNI_intervalCensoring"))
    }
    if (!all(c(y, rcen) > 0)) {
      stop("data should be positive")
    }
    if (!all(c(weight, rcenweight) >= 0)) {
      stop("weights should be non-negative")
    }
    if (r < 1 && any(methods == "RK")) {
      stop("sub-sampling available for UNI and PADE methods")
    }
    if (length(burn) != 1 || !is.finite(burn) || burn < 0 || burn > 100) {
      stop("burn should be a number between 0 and 100")
    }
    if (burn > 0 && has_censoring) {
      stop("burn is only available for uncensored data")
    }
    plot_enabled <- isTRUE(plot_progress)
    use_km_plot <- rightCensored && has_censoring
    if (plot_enabled) {
      if (!is.vector(y) || is.matrix(y) || !is.numeric(y)) {
        stop("plot_progress requires y to be a numeric vector")
      }
      if (length(plot_every) != 1 || !is.finite(plot_every) || plot_every < 1) {
        stop("plot_every should be a positive integer")
      }
      if (length(plot_breaks) != 1 || !is.finite(plot_breaks) || plot_breaks < 1) {
        stop("plot_breaks should be a positive integer")
      }
      if (length(plot_pause) != 1 || !is.finite(plot_pause) || plot_pause < 0) {
        stop("plot_pause should be a non-negative number")
      }
      plot_every <- as.integer(plot_every)
      plot_breaks <- as.integer(plot_breaks)
      plot_pause <- as.numeric(plot_pause)
      y_plot <- as.numeric(y)
      y_plot <- y_plot[is.finite(y_plot)]
      if (length(y_plot) < 2) {
        stop("plot_progress requires at least two finite observations")
      }
      if (use_km_plot) {
        rcen_plot <- as.numeric(rcen)
        rcen_plot <- rcen_plot[is.finite(rcen_plot)]
        w_plot <- if (length(weight) == 0) rep(1, length(y_plot)) else as.numeric(weight)
        rcw_plot <- if (length(rcenweight) == 0) rep(1, length(rcen_plot)) else as.numeric(rcenweight)
        if (length(w_plot) != length(y_plot) || length(rcw_plot) != length(rcen_plot)) {
          stop("incompatible plotting inputs for right-censored data")
        }
        km_raw <- data.frame(
          time = c(y_plot, rcen_plot),
          weight = c(w_plot, rcw_plot),
          event = c(rep(1, length(y_plot)), rep(0, length(rcen_plot)))
        )
        km_total <- stats::aggregate(weight ~ time, data = km_raw, FUN = sum)
        km_event <- stats::aggregate(weight ~ time, data = km_raw[km_raw$event == 1, , drop = FALSE], FUN = sum)
        event_weight <- numeric(nrow(km_total))
        if (nrow(km_event) > 0) {
          event_weight[match(km_event$time, km_total$time)] <- km_event$weight
        }
        risk <- sum(km_total$weight)
        surv <- numeric(nrow(km_total))
        s_now <- 1
        for (i in seq_len(nrow(km_total))) {
          if (risk > 0 && event_weight[i] > 0) {
            s_now <- s_now * (1 - event_weight[i] / risk)
          }
          surv[i] <- s_now
          risk <- risk - km_total$weight[i]
        }
        km_plot <- list(
          time = c(0, km_total$time),
          surv = c(1, pmax(0, pmin(1, surv)))
        )
        x_data <- c(y_plot, rcen_plot)
      } else {
        x_data <- y_plot
      }
      x_left <- min(x_data, na.rm = TRUE)
      x_right <- max(x_data, na.rm = TRUE)
      if (!is.finite(x_left) || !is.finite(x_right)) {
        stop("plot_progress requires finite plotting range")
      }
      if (use_km_plot) {
        x_left <- 0
      }
      if (x_right <= x_left) {
        pad <- max(abs(x_left) * 0.05, 1e-8)
        x_left <- x_left - pad
        x_right <- x_right + pad
      }
      grid_plot <- seq(x_left, x_right, length.out = 300)
      ll_start <- max(1L, as.integer(plot_every))
      plot_state <- function(model_obj, iter, ll_hist) {
        op <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(op), add = TRUE)
        graphics::par(mfrow = c(1, 2))
        if (use_km_plot) {
          graphics::plot(
            km_plot$time,
            km_plot$surv,
            type = "s",
            col = "black",
            xlim = c(x_left, x_right),
            ylim = c(0, 1),
            xlab = "",
            ylab = "",
            main = "",
            xaxs = "i"
          )
          surv_fit <- tryCatch(
            cdf(model_obj, grid_plot, lower.tail = FALSE),
            error = function(e) rep(NA_real_, length(grid_plot))
          )
          if (all(is.finite(surv_fit))) {
            graphics::lines(grid_plot, surv_fit, col = "red", lwd = 2)
          }
        } else {
          hobj <- graphics::hist(
            y_plot,
            breaks = plot_breaks,
            plot = FALSE
          )
          dfit <- tryCatch(
            dens(model_obj, grid_plot),
            error = function(e) rep(NA_real_, length(grid_plot))
          )
          y_top <- max(hobj$density, na.rm = TRUE)
          if (all(is.finite(dfit))) {
            y_top <- max(y_top, dfit, na.rm = TRUE)
          }
          if (!is.finite(y_top) || y_top <= 0) {
            y_top <- 1
          }
          y_top <- 1.05 * y_top
          graphics::hist(
            y_plot,
            probability = TRUE,
            breaks = plot_breaks,
            col = "grey90",
            border = "white",
            ylim = c(0, y_top),
            xlab = "",
            ylab = "",
            main = ""
          )
          if (all(is.finite(dfit))) {
            graphics::lines(grid_plot, dfit, col = "red", lwd = 2)
          }
        }
        if (iter < ll_start) {
          graphics::plot.new()
          graphics::box()
        } else {
          idx <- seq.int(ll_start, iter)
          graphics::plot(
            idx,
            ll_hist[idx],
            type = "l",
            col = "grey40",
            lwd = 2,
            xlab = "",
            ylab = "",
            main = ""
          )
        }
        utils::flush.console()
        if (plot_pause > 0) {
          Sys.sleep(plot_pause)
        }
      }
    }
    is_iph <- methods::is(x, "iph")
    merlang_fit <- .merlang_fit_blocks(x, merlang_blocks)
    erlang_fit <- if (is.null(erlang)) {
      grepl("\\berlang\\b", tolower(x@name))
    } else {
      isTRUE(erlang)
    }
    if (length(merlang_fit) > 0 && erlang_fit) {
      stop("use either erlang or merlang_blocks, not both")
    }
    if (is_iph) {
      par_g <- x@gfun$pars
      inv_g <- x@gfun$inverse
    }
    LL <- if (is_iph) {
      if(rightCensored){
        eval(parse(text = paste("logLikelihoodM", x@gfun$name, "_", methods[2], sep = "")))
      }else{
        eval(parse(text = paste("logLikelihoodM", x@gfun$name, "_UNI_intervalCensoring", sep = "")))
      }
    } else {
      if(rightCensored){
        eval(parse(text = paste("logLikelihoodPH_", methods[2], sep = "")))
      }else{
        eval(parse(text = "logLikelihoodPH_UNI_intervalCensoring"))
      }
    }
    
    A <- data_aggregation(y, weight)
    y_full <- A$un_obs
    weight_full <- A$weights
    
    if (length(rcen) > 0) {
      B <- data_aggregation(rcen, rcenweight)
      rcen_full <- B$un_obs
      rcenweight_full <- B$weights
    } else {
      rcen_full <- rcen
      rcenweight_full <- rcenweight
    }
    
    draw_subsample <- function(obs, obs_weight) {
      if (r >= 1) {
        return(list(obs = obs, weight = obs_weight))
      }
      n <- if (is.matrix(obs)) nrow(obs) else length(obs)
      if (n == 0) {
        return(list(obs = obs, weight = obs_weight))
      }
      size <- max(1L, floor(r * n))
      idx <- sample.int(n, size = size)
      if (is.matrix(obs)) {
        list(obs = obs[idx, , drop = FALSE], weight = obs_weight[idx])
      } else {
        list(obs = obs[idx], weight = obs_weight[idx])
      }
    }
    
    burn_steps <- floor(stepsEM * burn / 100)
    y_burn_full <- y_full
    weight_burn_full <- weight_full
    if (burn_steps > 0) {
      htmp <- graphics::hist(as.numeric(y_full), plot = FALSE)
      bin_id <- cut(
        as.numeric(y_full),
        breaks = htmp$breaks,
        include.lowest = TRUE,
        labels = FALSE
      )
      bin_w <- tapply(as.numeric(weight_full), bin_id, sum)
      keep <- which(!is.na(bin_w) & bin_w > 0)
      if (length(keep) > 0) {
        y_burn_full <- htmp$mids[keep]
        weight_burn_full <- as.numeric(bin_w[keep])
      } else {
        burn_steps <- 0
      }
    }
    
    get_iteration_data <- function(k) {
      if (burn_steps > 0 && k <= burn_steps) {
        y_pool <- y_burn_full
        weight_pool <- weight_burn_full
      } else {
        y_pool <- y_full
        weight_pool <- weight_full
      }
      obs_draw <- draw_subsample(y_pool, weight_pool)
      if (r < 1 && length(rcen_full) > 0) {
        rcen_draw <- draw_subsample(rcen_full, rcenweight_full)
        rcen_obs <- rcen_draw$obs
        rcen_w <- rcen_draw$weight
      } else {
        rcen_obs <- rcen_full
        rcen_w <- rcenweight_full
      }
      list(y = obs_draw$obs, weight = obs_draw$weight, rcen = rcen_obs, rcenweight = rcen_w)
    }
    
    ph_par <- x@pars
    alpha_fit <- clone_vector(ph_par$alpha)
    S_fit <- clone_matrix(ph_par$S)
    track <- rep(NA_real_, stepsEM)
    final_ll <- NA_real_
    
    options(digits.secs = 4)
    cat(format(Sys.time(), format = "%H:%M:%OS"), ": EM started", sep = "")
    cat("\n", sep = "")
    if (plot_enabled) {
      plot_state(x, iter = 0, ll_hist = track)
    }
    
    if (!is_iph) {
      for (k in seq_len(stepsEM)) {
        iter_data <- get_iteration_data(k)
        y <- iter_data$y
        weight <- iter_data$weight
        rcen <- iter_data$rcen
        rcenweight <- iter_data$rcenweight
        
        epsilon1 <- switch(which(methods[1] == c("RK", "UNI", "PADE")),
                           if (!is.na(rkstep)) rkstep else default_step_length(S_fit),
                           if (!is.na(uni_epsilon)) uni_epsilon else 1e-4,
                           0
        )
        epsilon2 <- switch(which(methods[2] == c("RK", "UNI", "PADE")),
                           if (!is.na(rkstep)) rkstep else default_step_length(S_fit),
                           if (!is.na(uni_epsilon)) uni_epsilon else 1e-4,
                           0
        )
        EMstep(epsilon1, alpha_fit, S_fit, y, weight, rcen, rcenweight, erlang_fit, merlang_fit)
        ll_cur <- LL(epsilon2, alpha_fit, S_fit, y, weight, rcen, rcenweight)
        track[k] <- ll_cur
        if (k %% every == 0) {
          cat("\r", "iteration:", k,
              ", logLik:", ll_cur,
              sep = " "
          )
        }
        if (plot_enabled && (k %% plot_every == 0)) {
          plot_state(ph(alpha = alpha_fit, S = S_fit), iter = k, ll_hist = track)
        }
      }
      x@pars$alpha <- alpha_fit
      x@pars$S <- S_fit
      if (length(merlang_fit) > 0) {
        S_attr <- x@pars$S
        attr(S_attr, "merlang_blocks") <- as.integer(merlang_fit)
        x@pars$S <- S_attr
      }
      epsilon_ll <- switch(
        which(methods[2] == c("RK", "UNI", "PADE")),
        if (!is.na(rkstep)) rkstep else default_step_length(S_fit),
        if (!is.na(uni_epsilon)) uni_epsilon else 1e-4,
        0
      )
      final_ll <- LL(epsilon_ll, alpha_fit, S_fit, y_full, weight_full, rcen_full, rcenweight_full)
      x@fit <- list(logLik = final_ll, nobs = sum(A$weights), logLikHist = track)
    }
    if (is_iph) {
      for (k in seq_len(stepsEM)) {
        iter_data <- get_iteration_data(k)
        y <- iter_data$y
        weight <- iter_data$weight
        rcen <- iter_data$rcen
        rcenweight <- iter_data$rcenweight
        trans_weight <- weight
        trans_rcenweight <- rcenweight
        
        if (x@gfun$name != "gev") {
          trans <- inv_g(par_g, y)
          trans_cens <- inv_g(par_g, rcen)
        } else {
          t <- inv_g(par_g, y, weight)
          tc <- inv_g(par_g, rcen, rcenweight)
          trans <- t$obs
          trans_weight <- t$weight
          trans_cens <- tc$obs
          trans_rcenweight <- tc$weight
        }
        epsilon1 <- switch(which(methods[1] == c("RK", "UNI", "PADE")),
                           if (!is.na(rkstep)) rkstep else default_step_length(S_fit),
                           if (!is.na(uni_epsilon)) uni_epsilon else 1e-4,
                           0
        )
        epsilon2 <- switch(which(methods[2] == c("RK", "UNI", "PADE")),
                           if (!is.na(rkstep)) rkstep else default_step_length(S_fit),
                           if (!is.na(uni_epsilon)) uni_epsilon else 1e-4,
                           0
        )
        EMstep(epsilon1, alpha_fit, S_fit, trans, trans_weight, trans_cens, trans_rcenweight, erlang_fit, merlang_fit)
        opt <- suppressWarnings(
          stats::optim(
            par = par_g,
            fn = LL,
            h = epsilon2,
            alpha = alpha_fit,
            S = S_fit,
            obs = y,
            weight = weight,
            rcens = rcen,
            rcweight = rcenweight,
            hessian = (k == stepsEM),
            method = optim_method,
            control = control
          )
        )
        track[k] <- opt$value
        par_g <- opt$par
        
        if (k %% every == 0) {
          cat("\r", "iteration:", k,
              ", logLik:", opt$value,
              sep = " "
          )
        }
        if (plot_enabled && (k %% plot_every == 0)) {
          plot_state(
            iph(ph(alpha = alpha_fit, S = S_fit), gfun = x@gfun$name, gfun_pars = par_g),
            iter = k,
            ll_hist = track
          )
        }
      }
      x@pars$alpha <- alpha_fit
      x@pars$S <- S_fit
      if (length(merlang_fit) > 0) {
        S_attr <- x@pars$S
        attr(S_attr, "merlang_blocks") <- as.integer(merlang_fit)
        x@pars$S <- S_attr
      }
      epsilon_ll <- switch(
        which(methods[2] == c("RK", "UNI", "PADE")),
        if (!is.na(rkstep)) rkstep else default_step_length(S_fit),
        if (!is.na(uni_epsilon)) uni_epsilon else 1e-4,
        0
      )
      final_ll <- LL(
        epsilon_ll,
        alpha_fit,
        S_fit,
        par_g,
        y_full,
        weight_full,
        rcen_full,
        rcenweight_full
      )
      x@fit <- list(logLik = final_ll, nobs = sum(A$weights), logLikHist = track)
      x <- iph(x, gfun = x@gfun$name, gfun_pars = par_g)
    }
    if (plot_enabled && (!is.finite(final_ll) || (stepsEM %% plot_every != 0))) {
      plot_state(x, iter = stepsEM, ll_hist = track)
    }
    
    cat("\n", format(Sys.time(), format = "%H:%M:%OS"), ": EM finalized", sep = "")
    cat("\n", sep = "")
    
    x
  }
)

#' Aggregate data to improve computation speed
#' 
#' @param y Observations. Either a vector or a matrix.
#' @param w Respective weights of observations
#' 
#' @return Returns a named list with unique observations and associated weights. If y is a vector then the unique observations are given, otherwise the unique rows are returned.
#' 
data_aggregation <- function(y, w) {
  if ((length(w) == 0) && (is.vector(y))) w <- rep(1, length(y))
  if ((length(w) == 0) && (is.matrix(y))) w <- rep(1, nrow(y))
  
  observations <- cbind(y, w)
  mat <- data.frame(observations)
  
  if(is.vector(y)){
    names(mat) <- c("obs", "weight")
    agg <- stats::aggregate(mat$weight,
                            by = list(un_obs = mat$obs),
                            FUN = sum)
  }else if(is.matrix(y)){
    names(mat) <- c("lower", "upper", "weight")
    agg <- stats::aggregate(mat$weight,
                            by = list(lower = mat$lower, upper = mat$upper),
                            FUN = sum)
    
    agg$un_obs <- as.matrix(agg[, c("lower", "upper")])
  }
  
  list(un_obs = agg$un_obs, weights = agg$x)
}

#' Loglikelihood method for ph class
#'
#' @param object An object of class \linkS4class{ph}.
#'
#' @return An object of class logLik.
#' @export
#'
#' @examples
#' obj <- iph(ph(structure = "general", dimension = 2), gfun = "weibull", gfun_pars = 2)
#' data <- sim(obj, n = 100)
#' fitted_ph <- fit(obj, data, stepsEM = 10)
#' logLik(fitted_ph)
setMethod("logLik", "ph", function(object) {
  ll <- object@fit$logLik
  attr(ll, "nobs") <- object@fit$nobs
  attr(ll, "df") <- sum(unlist(coef(object)) != 0) - 1
  class(ll) <- "logLik"
  ll
})

#' Coef method for ph class
#'
#' @param object An object of class \linkS4class{ph}.
#'
#' @return Parameters of ph model.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' coef(obj)
setMethod("coef", c(object = "ph"), function(object) {
  object@pars
})

#' LRT method for ph class
#'
#' @param x,y Objects of class \linkS4class{ph}.
#'
#' @return LRT between the models.
#' @export
#' @importFrom stats pchisq
#'
setMethod("LRT", c(x = "ph", y = "ph"), function(x, y) {
  LR <- 2 * abs(logLik(y) - logLik(x))
  degrees <- abs(attributes(logLik(y))$df - attributes(logLik(x))$df)
  c(LR = LR, p.val = pchisq(LR, df = degrees, lower.tail = FALSE))
})

#' TVR method for ph class
#'
#' @param x An object of class \linkS4class{ph}.
#' @param rew A vector of rewards.
#'
#' @return An object of the of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' obj <- ph(structure = "general")
#' TVR(obj, c(1, 2, 3))
setMethod("TVR", c(x = "ph"), function(x, rew) {
  if (length(x@pars$alpha) != length(rew)) {
    stop("vector of rewards of wrong dimension")
  }
  if (any(rew < 0)) {
    stop("vector of rewards with negative entries")
  }
  mar_par <- tvr_ph(x@pars$alpha, x@pars$S, rew)
  ph(alpha = mar_par[[1]], S = mar_par[[2]])
})

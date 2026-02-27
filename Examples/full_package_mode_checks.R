#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(matrixdist))
set.seed(20260227)

.passes <- character(0)
.fails <- character(0)
.skips <- character(0)

assert_true <- function(cond, msg = "assertion failed") {
  if (!isTRUE(cond)) {
    stop(msg, call. = FALSE)
  }
}

assert_finite <- function(x, msg = "non-finite values found") {
  if (any(!is.finite(as.numeric(x)))) {
    stop(msg, call. = FALSE)
  }
}

assert_prob <- function(alpha, tol = 1e-6) {
  assert_true(all(alpha >= -tol), "alpha has negative entries")
  assert_true(abs(sum(alpha) - 1) <= 1e-3, "alpha does not sum to 1")
}

safe_r <- function(S, frac = 0.5, min_r = 1e-4) {
  bounds <- if (is.list(S)) {
    vapply(S, function(M) -max(Re(eigen(M, only.values = TRUE)$values)), numeric(1))
  } else {
    -max(Re(eigen(S, only.values = TRUE)$values))
  }
  b <- min(bounds)
  if (!is.finite(b) || b <= min_r) {
    return(min_r)
  }
  frac * b
}

starts_from_blocks <- function(block_sizes) {
  cumsum(c(1L, block_sizes))[seq_along(block_sizes)]
}

run_check <- function(name, expr) {
  cat(sprintf("[CHECK] %s ... ", name))
  t0 <- proc.time()[3]
  ok <- TRUE
  err <- NULL
  tryCatch(
    force(expr),
    error = function(e) {
      ok <<- FALSE
      err <<- conditionMessage(e)
    }
  )
  elapsed <- round(proc.time()[3] - t0, 2)
  if (ok) {
    .passes <<- c(.passes, name)
    cat(sprintf("OK (%.2fs)\n", elapsed))
  } else {
    .fails <<- c(.fails, sprintf("%s :: %s", name, err))
    cat("FAIL\n")
  }
}

skip_check <- function(name, reason) {
  .skips <<- c(.skips, sprintf("%s :: %s", name, reason))
  cat(sprintf("[CHECK] %s ... SKIP (%s)\n", name, reason))
}

run_check("ph: structure constructors + base methods + fit(UNI)", {
  structures <- c("general", "hyperexponential", "coxian", "gcoxian", "gerlang", "erlang", "erland")
  for (st in structures) {
    x <- ph(structure = st, dimension = 4)
    y <- sim(x, n = 250)
    assert_true(all(y > 0), paste("non-positive simulation for", st))
    q <- as.numeric(stats::quantile(y, probs = c(0.2, 0.5, 0.8)))
    assert_finite(dens(x, q), paste("dens failed for", st))
    assert_finite(cdf(x, q), paste("cdf failed for", st))
    assert_finite(haz(x, q), paste("haz failed for", st))
    assert_finite(quan(x, c(0.25, 0.5, 0.75)), paste("quan failed for", st))

    w <- stats::runif(length(y), 0.9, 1.1)
    f <- fit(
      ph(structure = st, dimension = 4),
      y = y,
      weight = w,
      stepsEM = 8,
      every = 4,
      methods = c("UNI", "UNI"),
      erlang = if (st %in% c("erlang", "erland")) TRUE else NULL
    )
    assert_prob(f@pars$alpha)
    assert_finite(f@pars$S, paste("non-finite S for", st))
    assert_finite(logLik(f), paste("non-finite logLik for", st))
  }
})

run_check("ph: RK/UNI/PADE + weights + right/interval censoring + subsampling", {
  truth <- ph(structure = "gerlang", dimension = 5)
  y <- sim(truth, n = 700)

  c_lim <- as.numeric(stats::quantile(y, 0.85))
  unc <- y <= c_lim
  y_unc <- y[unc]
  y_rc <- rep(c_lim, sum(!unc))
  w_unc <- stats::runif(length(y_unc), 0.8, 1.2)
  w_rc <- stats::runif(length(y_rc), 0.8, 1.2)

  f_uni <- fit(
    ph(structure = "gerlang", dimension = 5),
    y = y_unc,
    weight = w_unc,
    rcen = y_rc,
    rcenweight = w_rc,
    stepsEM = 12,
    every = 6,
    methods = c("UNI", "UNI"),
    r = 0.75
  )
  assert_finite(logLik(f_uni), "UNI right-censored fit failed")

  f_pade <- fit(
    ph(structure = "gerlang", dimension = 5),
    y = y_unc,
    weight = w_unc,
    rcen = y_rc,
    rcenweight = w_rc,
    stepsEM = 12,
    every = 6,
    methods = c("PADE", "PADE"),
    r = 0.75
  )
  assert_finite(logLik(f_pade), "PADE right-censored fit failed")

  f_rk <- fit(
    ph(structure = "gerlang", dimension = 5),
    y = y_unc,
    weight = w_unc,
    rcen = y_rc,
    rcenweight = w_rc,
    stepsEM = 12,
    every = 6,
    methods = c("RK", "RK")
  )
  assert_finite(logLik(f_rk), "RK right-censored fit failed")

  idx_int <- sample(seq_along(y), size = 140)
  width <- stats::runif(length(idx_int), 0.05, 0.25)
  rc_int <- cbind(pmax(1e-8, y[idx_int] - width), y[idx_int] + width)
  y_obs <- y[-idx_int]
  w_obs <- stats::runif(length(y_obs), 0.9, 1.1)
  w_int <- stats::runif(nrow(rc_int), 0.9, 1.1)

  f_int <- fit(
    ph(structure = "coxian", dimension = 5),
    y = y_obs,
    weight = w_obs,
    rcen = rc_int,
    rcenweight = w_int,
    stepsEM = 10,
    every = 5,
    methods = c("UNI", "UNI")
  )
  assert_finite(logLik(f_int), "interval-censored UNI fit failed")
})

run_check("ph: exact Erlang constrained M-step", {
  p <- 5
  lambda_true <- 0.7
  alpha <- c(1, rep(0, p - 1))
  S <- matrix(0, p, p)
  diag(S) <- -lambda_true
  S[cbind(1:(p - 1), 2:p)] <- lambda_true
  truth <- ph(alpha = alpha, S = S)

  y <- sim(truth, n = 600)
  fit_e <- fit(
    ph(structure = "erlang", dimension = p),
    y = y,
    stepsEM = 40,
    every = 20,
    methods = c("UNI", "UNI"),
    erlang = TRUE
  )

  sup <- fit_e@pars$S[cbind(1:(p - 1), 2:p)]
  assert_true(stats::sd(sup) < 1e-9, "superdiagonal should be constant for Erlang")
  assert_true(abs(fit_e@pars$alpha[1] - 1) < 1e-12, "Erlang alpha[1] should be 1")
  assert_true(max(abs(fit_e@pars$alpha[-1])) < 1e-12, "Erlang alpha[-1] should be 0")
})

run_check("ph: merlang constructors + UNI/RK/PADE exact block updates", {
  blocks <- c(2, 3, 2)
  starts <- starts_from_blocks(blocks)
  via_ph <- ph(structure = "merlang", block_sizes = blocks)
  assert_true(methods::is(via_ph, "ph"), "ph(structure='merlang') should return ph")
  via_ph_alias <- ph(structure = "mixederlang", block_sizes = blocks)
  assert_true(methods::is(via_ph_alias, "ph"), "ph(structure='mixederlang') should return ph")
  truth <- merlang_ph(block_sizes = blocks, probs = c(0.25, 0.5, 0.25), rates = c(0.6, 1.1, 1.8))
  y <- sim(truth, n = 800)

  fit_uni <- fit(
    merlang_ph(block_sizes = blocks),
    y = y,
    stepsEM = 50,
    every = 10,
    methods = c("UNI", "UNI"),
    merlang_blocks = blocks
  )
  fit_rk <- fit(
    merlang_ph(block_sizes = blocks),
    y = y,
    stepsEM = 20,
    every = 10,
    methods = c("RK", "RK"),
    merlang_blocks = blocks
  )
  fit_pade <- fit(
    merlang_ph(block_sizes = blocks),
    y = y,
    stepsEM = 20,
    every = 10,
    methods = c("PADE", "PADE"),
    merlang_blocks = blocks
  )

  for (obj in list(fit_uni, fit_rk, fit_pade)) {
    diag_rates <- -diag(obj@pars$S)[starts]
    fwd_rates <- obj@pars$S[cbind(starts, starts + 1L)]
    assert_true(max(abs(diag_rates - fwd_rates)) < 1e-9, "block forward rates must match block diagonal rates")
    assert_true(identical(as.integer(attr(obj@pars$S, "merlang_blocks")), as.integer(blocks)), "merlang_blocks attribute not retained")
    assert_prob(obj@pars$alpha)
  }
})

run_check("iph: multiple inhomogeneity families (weibull/pareto/lognormal/gompertz)", {
  alpha0 <- c(1, 0, 0, 0)
  S0 <- matrix(0, 4, 4)
  diag(S0) <- -c(1.2, 1.0, 0.9, 0.8)
  S0[cbind(1:3, 2:4)] <- c(0.7, 0.6, 0.5)
  gfuns <- list(
    list(name = "weibull", pars = 1.4),
    list(name = "pareto", pars = 1.2),
    list(name = "lognormal", pars = 1.1),
    list(name = "gompertz", pars = 0.02)
  )

  for (g in gfuns) {
    x <- iph(ph(alpha = alpha0, S = S0), gfun = g$name, gfun_pars = g$pars)
    y <- sim(x, n = 320)
    w <- stats::runif(length(y), 0.9, 1.1)
    f <- fit(
      iph(ph(alpha = alpha0, S = S0), gfun = g$name, gfun_pars = g$pars),
      y = y,
      weight = w,
      stepsEM = 8,
      every = 4,
      methods = c("UNI", "UNI")
    )
    assert_finite(logLik(f), paste("iph fit failed for", g$name))
    assert_finite(f@pars$S, paste("iph S invalid for", g$name))
  }
})

run_check("iph: exact Erlang + exact merlang in fit()", {
  p <- 5
  x_e <- iph(ph(structure = "erlang", dimension = p), gfun = "weibull", gfun_pars = 1.2)
  y_e <- sim(x_e, n = 500)
  f_e <- fit(
    iph(ph(structure = "erlang", dimension = p), gfun = "weibull", gfun_pars = 1.0),
    y = y_e,
    stepsEM = 20,
    every = 10,
    methods = c("UNI", "UNI"),
    erlang = TRUE
  )
  sup <- f_e@pars$S[cbind(1:(p - 1), 2:p)]
  assert_true(stats::sd(sup) < 1e-9, "inhomogeneous Erlang fit is not exactly constrained")

  blocks <- c(2, 3)
  starts <- starts_from_blocks(blocks)
  x_m <- iph(merlang_ph(block_sizes = blocks), gfun = "weibull", gfun_pars = 1.1)
  y_m <- sim(x_m, n = 550)
  f_m <- fit(
    iph(merlang_ph(block_sizes = blocks), gfun = "weibull", gfun_pars = 1.0),
    y = y_m,
    stepsEM = 25,
    every = 5,
    methods = c("UNI", "UNI"),
    merlang_blocks = blocks
  )
  assert_true(identical(as.integer(attr(f_m@pars$S, "merlang_blocks")), as.integer(blocks)), "merlang attr missing after iph fit")
  assert_true(max(abs((-diag(f_m@pars$S)[starts]) - f_m@pars$S[cbind(starts, starts + 1L)])) < 1e-9, "inhomogeneous merlang fit not block-constrained")
})

run_check("reg + inhomogeneity_reg: erlang/merlang, censoring, weights, sph tools", {
  p <- 5
  x_true <- iph(ph(structure = "erlang", dimension = p), gfun = "weibull", gfun_pars = 1.25)
  y_full <- sim(x_true, n = 700)

  c_lim <- as.numeric(stats::quantile(y_full, 0.85))
  unc <- y_full <= c_lim
  y <- y_full[unc]
  rc <- rep(c_lim, sum(!unc))
  X_full <- matrix(stats::rnorm(length(y_full)), ncol = 1)
  X <- rbind(X_full[unc, , drop = FALSE], X_full[!unc, , drop = FALSE])
  w <- stats::runif(length(y), 0.9, 1.1)
  w_rc <- stats::runif(length(rc), 0.9, 1.1)

  fit_reg_e <- reg(
    x = iph(ph(structure = "erlang", dimension = p), gfun = "weibull", gfun_pars = 1.0),
    y = y,
    weight = w,
    rcen = rc,
    rcenweight = w_rc,
    X = X,
    stepsEM = 15,
    every = 5,
    methods = c("UNI", "UNI"),
    erlang = TRUE
  )
  sup <- fit_reg_e@pars$S[cbind(1:(p - 1), 2:p)]
  assert_true(stats::sd(sup) < 1e-9, "reg Erlang fit not exactly constrained")

  blocks <- c(2, 3)
  starts <- starts_from_blocks(blocks)
  x_true_m <- iph(merlang_ph(block_sizes = blocks), gfun = "weibull", gfun_pars = 1.2)
  y_m <- sim(x_true_m, n = 750)
  X_m <- matrix(stats::rnorm(length(y_m)), ncol = 1)

  fit_reg_m <- reg(
    x = iph(merlang_ph(block_sizes = blocks), gfun = "weibull", gfun_pars = 1.0),
    y = y_m,
    X = X_m,
    stepsEM = 20,
    every = 10,
    methods = c("UNI", "UNI"),
    merlang_blocks = blocks
  )
  assert_true(identical(as.integer(attr(fit_reg_m@pars$S, "merlang_blocks")), as.integer(blocks)), "reg merlang attr missing")
  assert_true(max(abs((-diag(fit_reg_m@pars$S)[starts]) - fit_reg_m@pars$S[cbind(starts, starts + 1L)])) < 1e-9, "reg merlang block constraints broken")

  prop_f <- function(theta, data) exp(drop(data %*% theta))
  inhom_f <- function(theta, data) exp(drop(data %*% theta))
  fit_inhom_m <- inhomogeneity_reg(
    x = iph(merlang_ph(block_sizes = blocks), gfun = "weibull", gfun_pars = 1.0),
    y = y_m,
    X = X_m,
    X2 = X_m,
    prop_f = prop_f,
    inhom_f = inhom_f,
    B0 = rep(0, 2),
    stepsEM = 12,
    every = 4,
    methods = c("UNI", "UNI"),
    merlang_blocks = blocks
  )
  assert_true(identical(as.integer(attr(fit_inhom_m@pars$S, "merlang_blocks")), as.integer(blocks)), "inhomogeneity_reg merlang attr missing")

  subj <- rep(0, length(fit_reg_e@coefs$B))
  ev <- evaluate(fit_reg_e, subject = subj)
  assert_true(methods::is(ev, "ph"), "evaluate(sph) should return ph")

  n_fisher <- min(80, length(y))
  fi <- Fisher(fit_reg_e, y = y[seq_len(n_fisher)], X = X[seq_len(n_fisher), , drop = FALSE], w = w[seq_len(n_fisher)])
  assert_finite(fi, "Fisher information contains non-finite values")
})

run_check("dph: structure constructors + methods + fit", {
  structures <- c("general", "hyperexponential", "coxian", "gcoxian", "gerlang", "erlang", "erland")
  for (st in structures) {
    x <- dph(structure = st, dimension = 4)
    y <- sim(x, n = 320)
    assert_true(all(y > 0), paste("non-positive dph simulation for", st))
    pts <- pmax(1L, as.integer(stats::quantile(y, probs = c(0.2, 0.5, 0.8))))
    assert_finite(dens(x, pts), paste("dph dens failed for", st))
    assert_finite(cdf(x, pts), paste("dph cdf failed for", st))
    assert_finite(pgf(x, c(0.3, 0.6)), paste("dph pgf failed for", st))

    w <- stats::runif(length(y), 0.9, 1.1)
    f <- fit(
      dph(structure = st, dimension = 4),
      y = y,
      weight = w,
      stepsEM = 12,
      every = 6,
      erlang = if (st %in% c("erlang", "erland")) TRUE else NULL
    )
    assert_prob(f@pars$alpha)
    assert_finite(f@pars$S, paste("dph fit produced invalid S for", st))
  }
})

run_check("dph: exact merlang blocks", {
  blocks <- c(3, 2)
  starts <- starts_from_blocks(blocks)
  via_dph <- dph(structure = "merlang", block_sizes = blocks)
  assert_true(methods::is(via_dph, "dph"), "dph(structure='merlang') should return dph")
  via_dph_alias <- dph(structure = "mixed_erlang", block_sizes = blocks)
  assert_true(methods::is(via_dph_alias, "dph"), "dph(structure='mixed_erlang') should return dph")
  truth <- merlang_dph(block_sizes = blocks, probs = c(0.45, 0.55), q = c(0.7, 0.85))
  y <- sim(truth, n = 900)
  w <- stats::runif(length(y), 0.95, 1.05)
  f <- fit(
    merlang_dph(block_sizes = blocks),
    y = y,
    weight = w,
    stepsEM = 60,
    every = 20,
    merlang_blocks = blocks
  )
  q_diag <- 1 - diag(f@pars$S)[starts]
  q_fwd <- f@pars$S[cbind(starts, starts + 1L)]
  assert_true(max(abs(q_diag - q_fwd)) < 1e-9, "merlang_dph block constraints broken")
  assert_true(identical(as.integer(attr(f@pars$S, "merlang_blocks")), as.integer(blocks)), "merlang_dph attr missing")
})

run_check("mph + miph: methods + fit (with censoring/subsampling)", {
  x <- mph(structure = c("general", "coxian"), dimension = 3)
  r0 <- safe_r(x@pars$S)
  eval_grid <- matrix(c(0.6 * r0, 0.5 * r0), ncol = 2)
  y <- sim(x, n = 260)
  assert_finite(dens(x, y[1:40, , drop = FALSE]), "mph dens failed")
  assert_finite(cdf(x, y[1:40, , drop = FALSE]), "mph cdf failed")
  assert_finite(laplace(x, eval_grid), "mph laplace failed")
  assert_finite(mgf(x, 0.25 * eval_grid), "mph mgf failed")
  assert_finite(moment(x, c(1, 1)), "mph moment failed")
  assert_finite(mean(x), "mph mean failed")
  assert_finite(var(x), "mph var failed")
  assert_finite(cor(x), "mph cor failed")

  delta <- matrix(1, nrow(y), ncol(y))
  delta[sample(length(delta), size = floor(0.12 * length(delta)))] <- 0
  fit_x <- fit(
    x,
    y = y,
    delta = delta,
    stepsEM = 8,
    r = r0
  )
  assert_prob(fit_x@pars$alpha)

  x_m <- miph(
    mph(structure = c("general", "coxian"), dimension = 3),
    gfun = c("weibull", "pareto"),
    gfun_pars = list(1.2, 1.1)
  )
  y_m <- sim(x_m, n = 220)
  x_m0 <- miph(
    mph(structure = c("general", "coxian"), dimension = 3),
    gfun = c("weibull", "pareto"),
    gfun_pars = list(1.1, 1.0)
  )
  fit_m <- fit(
    x_m0,
    y = y_m,
    stepsEM = 6,
    maxit = 20,
    r = safe_r(x_m0@pars$S)
  )
  assert_finite(unlist(fit_m@gfun$pars), "miph gfun parameters invalid after fit")
})

run_check("bivph + biviph: methods + fit", {
  x <- bivph(dimensions = c(3, 3))
  y <- sim(x, n = 240)
  assert_finite(dens(x, y[1:40, , drop = FALSE]), "bivph dens failed")
  assert_finite(moment(x, c(1, 1)), "bivph moment failed")
  assert_finite(mean(x), "bivph mean failed")
  assert_finite(var(x), "bivph var failed")
  assert_finite(cor(x), "bivph cor failed")
  assert_finite(laplace(x, matrix(c(0.7, 0.9), ncol = 2)), "bivph laplace failed")
  assert_finite(mgf(x, matrix(c(0.2, 0.1), ncol = 2)), "bivph mgf failed")
  assert_true(methods::is(marginal(x, 1), "ph"), "bivph marginal should be ph")
  assert_true(methods::is(linCom(x, c(1, 1)), "ph"), "bivph linCom should be ph")

  w <- stats::runif(nrow(y), 0.9, 1.1)
  f <- fit(bivph(dimensions = c(3, 3)), y = y, weight = w, stepsEM = 10, every = 5)
  assert_finite(f@pars$S11, "bivph fit invalid S11")

  x_i <- biviph(
    bivph(dimensions = c(3, 3)),
    gfun = c("weibull", "pareto"),
    gfun_pars = list(1.2, 1.1)
  )
  y_i <- sim(x_i, n = 220)
  f_i <- fit(
    biviph(
      bivph(dimensions = c(3, 3)),
      gfun = c("weibull", "pareto"),
      gfun_pars = list(1.1, 1.0)
    ),
    y = y_i,
    stepsEM = 6,
    every = 3
  )
  assert_finite(unlist(f_i@gfun$pars), "biviph gfun parameters invalid")
})

run_check("mdph + bivdph: methods + fit", {
  x_md <- mdph(structure = c("general", "coxian"), dimension = 3)
  y_md <- sim(x_md, n = 240)
  assert_finite(dens(x_md, y_md[1:40, , drop = FALSE]), "mdph dens failed")
  assert_finite(pgf(x_md, matrix(c(0.5, 0.4), ncol = 2)), "mdph pgf failed")
  assert_finite(moment(x_md, c(1, 1)), "mdph moment failed")
  assert_finite(mean(x_md), "mdph mean failed")
  assert_finite(var(x_md), "mdph var failed")
  assert_finite(cor(x_md), "mdph cor failed")
  f_md <- fit(x_md, y = y_md, weight = stats::runif(nrow(y_md), 0.9, 1.1), stepsEM = 10, every = 5)
  assert_finite(f_md@pars$S[[1]], "mdph fit invalid")

  x_bd <- bivdph(dimensions = c(3, 3))
  y_bd <- sim(x_bd, n = 240)
  assert_finite(dens(x_bd, y_bd[1:40, , drop = FALSE]), "bivdph dens failed")
  assert_finite(moment(x_bd, c(1, 1)), "bivdph moment failed")
  assert_finite(mean(x_bd), "bivdph mean failed")
  assert_finite(var(x_bd), "bivdph var failed")
  assert_finite(cor(x_bd), "bivdph cor failed")
  assert_finite(pgf(x_bd, matrix(c(0.5, 0.6), ncol = 2)), "bivdph pgf failed")
  f_bd <- fit(x_bd, y = y_bd, weight = stats::runif(nrow(y_bd), 0.9, 1.1), stepsEM = 10, every = 5)
  assert_finite(f_bd@pars$S11, "bivdph fit invalid S11")
})

run_check("MPHstar: methods + fit", {
  x <- MPHstar(structure = "general", dimension = 4, variables = 3)
  y <- sim(x, n = 260)
  assert_finite(mean(x), "MPHstar mean failed")
  assert_finite(var(x), "MPHstar var failed")
  assert_finite(cor(x), "MPHstar cor failed")
  assert_true(methods::is(marginal(x, 1), "ph"), "MPHstar marginal should be ph")
  assert_true(methods::is(linCom(x, c(1, 0, 1)), "ph"), "MPHstar linCom should be ph")

  w <- stats::runif(nrow(y), 0.9, 1.1)
  f <- fit(
    MPHstar(structure = "general", dimension = 4, variables = 3),
    y = y,
    weight = matrix(w, nrow = nrow(y), ncol = ncol(y)),
    stepsEM = 8,
    every = 4,
    r = 0.8
  )
  assert_finite(f@pars$S, "MPHstar fitted S invalid")
  assert_true(all(f@pars$R >= 0), "MPHstar rewards should remain non-negative")
})

if (requireNamespace("nnet", quietly = TRUE) && requireNamespace("reshape2", quietly = TRUE)) {
  run_check("MoE: ph, dph, mdph, bivdph", {
    n <- 80

    x_ph <- iph(ph(structure = "gerlang", dimension = 3), gfun = "weibull", gfun_pars = 1.2)
    y_ph <- sim(x_ph, n = n)
    data_ph <- data.frame(y = y_ph, x1 = stats::runif(n), x2 = stats::rnorm(n))
    dlt <- stats::rbinom(n, 1, 0.8)
    moe_ph <- MoE(x_ph, formula = y ~ x1 + x2, data = data_ph, delta = dlt, stepsEM = 2, every = 1)
    assert_true(is.list(moe_ph), "MoE(ph) should return a list")

    x_d <- dph(structure = "gerlang", dimension = 3)
    y_d <- sim(x_d, n = n)
    data_d <- data.frame(y = y_d, x1 = stats::runif(n), x2 = stats::rnorm(n))
    moe_d <- MoE(x_d, formula = y ~ x1 + x2, data = data_d, stepsEM = 2, every = 1)
    assert_true(is.list(moe_d), "MoE(dph) should return a list")

    x_md <- mdph(structure = c("general", "coxian"), dimension = 3)
    y_md <- sim(x_md, n = n)
    data_md <- data.frame(y1 = y_md[, 1], x1 = stats::runif(n), x2 = stats::rnorm(n))
    moe_md <- MoE(x_md, formula = y1 ~ x1 + x2, y = y_md, data = data_md, stepsEM = 2, every = 1)
    assert_true(is.list(moe_md), "MoE(mdph) should return a list")

    x_bd <- bivdph(dimensions = c(3, 3))
    y_bd <- sim(x_bd, n = n)
    data_bd <- data.frame(y1 = y_bd[, 1], x1 = stats::runif(n), x2 = stats::rnorm(n))
    moe_bd <- MoE(x_bd, formula = y1 ~ x1 + x2, y = y_bd, data = data_bd, stepsEM = 2, every = 1)
    assert_true(is.list(moe_bd), "MoE(bivdph) should return a list")
  })
} else {
  skip_check("MoE: ph, dph, mdph, bivdph", "nnet/reshape2 not available")
}

cat("\n============================================================\n")
cat(sprintf("Checks passed: %d\n", length(.passes)))
cat(sprintf("Checks failed: %d\n", length(.fails)))
cat(sprintf("Checks skipped: %d\n", length(.skips)))

if (length(.skips) > 0) {
  cat("\nSkipped checks:\n")
  for (s in .skips) cat(" - ", s, "\n", sep = "")
}

if (length(.fails) > 0) {
  cat("\nFailed checks:\n")
  for (f in .fails) cat(" - ", f, "\n", sep = "")
  stop("At least one package check failed.", call. = FALSE)
}

cat("\nAll requested package mode checks completed successfully.\n")

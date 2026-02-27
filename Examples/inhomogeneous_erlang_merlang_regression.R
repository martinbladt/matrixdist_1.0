#!/usr/bin/env Rscript

# Inhomogeneous PH regression examples with exact Erlang / mixture-of-Erlangs EM.

library(matrixdist)
set.seed(20260227)

extract_block_rates <- function(S, block_sizes) {
  starts <- cumsum(c(1L, block_sizes))[seq_along(block_sizes)]
  list(
    starts = starts,
    diag = -diag(S)[starts],
    forward = S[cbind(starts, starts + 1L)]
  )
}

cat("=== Inhomogeneous Erlang regression (reg) ===\n")
p <- 5
x_true_erlang <- iph(
  ph(structure = "erlang", dimension = p),
  gfun = "weibull",
  gfun_pars = 1.3
)
y_erlang <- sim(x_true_erlang, n = 700)
X_erlang <- matrix(stats::rnorm(length(y_erlang)), ncol = 1)

fit_erlang_reg <- reg(
  x = iph(ph(structure = "erlang", dimension = p), gfun = "weibull", gfun_pars = 1.0),
  y = y_erlang,
  X = X_erlang,
  stepsEM = 40,
  every = 10,
  methods = c("UNI", "UNI"),
  erlang = TRUE
)

sup_erlang <- fit_erlang_reg@pars$S[cbind(1:(p - 1), 2:p)]
cat("  superdiag values:", paste(round(sup_erlang, 5), collapse = ", "), "\n")
cat("  superdiag sd:", round(stats::sd(sup_erlang), 8), "\n")

cat("\n=== Inhomogeneous mixture-of-Erlangs regression (reg) ===\n")
block_sizes <- c(2, 3)
x_true_merlang <- iph(
  merlang_ph(block_sizes = block_sizes, probs = c(0.4, 0.6), rates = c(0.75, 1.35)),
  gfun = "weibull",
  gfun_pars = 1.2
)
y_merlang <- sim(x_true_merlang, n = 900)
X_merlang <- matrix(stats::rnorm(length(y_merlang)), ncol = 1)

fit_merlang_reg <- reg(
  x = iph(merlang_ph(block_sizes = block_sizes), gfun = "weibull", gfun_pars = 1.0),
  y = y_merlang,
  X = X_merlang,
  stepsEM = 50,
  every = 10,
  methods = c("UNI", "UNI"),
  merlang_blocks = block_sizes
)

rates_reg <- extract_block_rates(fit_merlang_reg@pars$S, block_sizes)
cat("  block weights:", paste(round(fit_merlang_reg@pars$alpha[rates_reg$starts], 4), collapse = ", "), "\n")
cat("  block rates:", paste(round(rates_reg$diag, 4), collapse = ", "), "\n")
cat("  block forward rates:", paste(round(rates_reg$forward, 4), collapse = ", "), "\n")

cat("\n=== Inhomogeneous mixture-of-Erlangs inhomogeneity_reg ===\n")
prop_f <- function(theta, data) exp(drop(data %*% theta))
inhom_f <- function(theta, data) exp(drop(data %*% theta))

fit_merlang_inhom <- inhomogeneity_reg(
  x = iph(merlang_ph(block_sizes = block_sizes), gfun = "weibull", gfun_pars = 1.0),
  y = y_merlang,
  X = X_merlang,
  X2 = X_merlang,
  prop_f = prop_f,
  inhom_f = inhom_f,
  B0 = rep(0, 2),
  stepsEM = 30,
  every = 10,
  methods = c("UNI", "UNI"),
  merlang_blocks = block_sizes
)

rates_inhom <- extract_block_rates(fit_merlang_inhom@pars$S, block_sizes)
cat("  block weights:", paste(round(fit_merlang_inhom@pars$alpha[rates_inhom$starts], 4), collapse = ", "), "\n")
cat("  block rates:", paste(round(rates_inhom$diag, 4), collapse = ", "), "\n")
cat("  block forward rates:", paste(round(rates_inhom$forward, 4), collapse = ", "), "\n")

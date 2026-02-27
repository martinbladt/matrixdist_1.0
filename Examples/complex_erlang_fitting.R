#!/usr/bin/env Rscript

# More complex examples for Erlang fitting:
# - weighted data
# - right-censoring
# - interval-censoring
# - larger number of EM iterations

library(matrixdist)
set.seed(20260227)

extract_lambda <- function(S) {
  p <- nrow(S)
  super <- S[cbind(1:(p - 1), 2:p)]
  list(super_mean = mean(super), diag_mean = mean(diag(S)))
}

extract_q <- function(S) {
  p <- nrow(S)
  super <- S[cbind(1:(p - 1), 2:p)]
  list(super_mean = mean(super), diag_mean = mean(diag(S)))
}

cat("=== Continuous PH (Erlang) with right-censoring + weights ===\n")
p <- 5
lambda_true <- 0.7
alpha_true <- c(1, rep(0, p - 1))
S_true <- matrix(0, p, p)
diag(S_true) <- -lambda_true
S_true[cbind(1:(p - 1), 2:p)] <- lambda_true
ph_true <- ph(alpha = alpha_true, S = S_true)

n <- 600
x <- sim(ph_true, n = n)

censor_limit <- as.numeric(stats::quantile(x, 0.90))
delta <- as.integer(x <= censor_limit)
y_unc <- x[delta == 1]
y_rc <- rep(censor_limit, sum(delta == 0))
w_unc <- stats::runif(length(y_unc), 0.8, 1.2)
w_rc <- stats::runif(length(y_rc), 0.8, 1.2)

ph_start <- ph(structure = "erlang", dimension = p)
ph_fit_rc <- fit(
  ph_start,
  y = y_unc,
  weight = w_unc,
  rcen = y_rc,
  rcenweight = w_rc,
  stepsEM = 80,
  every = 20,
  methods = c("UNI", "UNI"),
  uni_epsilon = 1e-4
)

ph_fit_info <- extract_lambda(ph_fit_rc@pars$S)
cat("true lambda:", lambda_true, "\n")
cat("fitted superdiag mean:", round(ph_fit_info$super_mean, 6), "\n")
cat("fitted diagonal mean:", round(ph_fit_info$diag_mean, 6), "\n")
cat("fitted S:\n")
print(round(ph_fit_rc@pars$S, 4))

cat("\n=== Continuous PH (Erlang) with interval-censoring ===\n")
idx_int <- sample(seq_len(n), size = floor(0.30 * n))
width <- stats::runif(length(idx_int), 0.05, 0.25)
lower <- pmax(1e-8, x[idx_int] - width)
upper <- x[idx_int] + width
rc_int <- cbind(lower, upper)

y_obs <- x[-idx_int]
w_obs <- rep(1, length(y_obs))
w_int <- rep(1, nrow(rc_int))

ph_fit_int <- fit(
  ph(structure = "erlang", dimension = p),
  y = y_obs,
  weight = w_obs,
  rcen = rc_int,
  rcenweight = w_int,
  stepsEM = 60,
  every = 15,
  methods = c("UNI", "UNI"),
  uni_epsilon = 1e-5
)

ph_int_info <- extract_lambda(ph_fit_int@pars$S)
cat("interval-censored fitted superdiag mean:", round(ph_int_info$super_mean, 6), "\n")
cat("interval-censored fitted diagonal mean:", round(ph_int_info$diag_mean, 6), "\n")

cat("\n=== Discrete DPH (Erlang) with weights and more EM iterations ===\n")
q_true <- 0.82
S_d_true <- matrix(0, p, p)
diag(S_d_true) <- 1 - q_true
S_d_true[cbind(1:(p - 1), 2:p)] <- q_true
dph_true <- dph(alpha = alpha_true, S = S_d_true)

z <- sim(dph_true, n = 800)
wz <- stats::runif(length(z), 0.9, 1.1)

dph_fit <- fit(
  dph(structure = "erlang", dimension = p),
  y = z,
  weight = wz,
  stepsEM = 120,
  every = 30
)

dph_fit_info <- extract_q(dph_fit@pars$S)
cat("true q:", q_true, "\n")
cat("fitted superdiag mean:", round(dph_fit_info$super_mean, 6), "\n")
cat("fitted diagonal mean:", round(dph_fit_info$diag_mean, 6), "\n")
cat("fitted S:\n")
print(round(dph_fit@pars$S, 4))

#!/usr/bin/env Rscript

# Minimal examples for the new Erlang features.

library(matrixdist)
set.seed(123)

cat("=== Initialization ===\n")
ph_erlang <- ph(structure = "erlang", dimension = 4)
ph_erland <- ph(structure = "erland", dimension = 4) # alias
cat("alpha (erlang):\n")
print(ph_erlang@pars$alpha)
cat("S (erlang):\n")
print(round(ph_erlang@pars$S, 4))
cat("S (erland alias):\n")
print(round(ph_erland@pars$S, 4))

cat("\n=== Continuous PH fit (auto Erlang detection) ===\n")
y <- sim(ph_erlang, n = 80)
ph_fit <- fit(ph_erlang, y, stepsEM = 3, every = 1)
print(round(ph_fit@pars$S, 4))

cat("\n=== Continuous PH fit (explicit erlang = TRUE) ===\n")
ph_custom <- ph(alpha = ph_erlang@pars$alpha, S = ph_erlang@pars$S)
ph_fit_forced <- fit(ph_custom, y, stepsEM = 3, every = 1, erlang = TRUE)
print(round(ph_fit_forced@pars$S, 4))

cat("\n=== Discrete DPH fit (auto Erlang detection) ===\n")
dph_erlang <- dph(structure = "erlang", dimension = 4)
yd <- sim(dph_erlang, n = 120)
dph_fit <- fit(dph_erlang, yd, stepsEM = 3, every = 1)
print(round(dph_fit@pars$S, 4))


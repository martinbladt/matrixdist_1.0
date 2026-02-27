#!/usr/bin/env Rscript

# Complex mixture-of-Erlangs fitting examples:
# - exact blockwise EM (continuous + discrete)
# - UNI, RK, and PADE for continuous PH
# - weighted and right-censored data

library(matrixdist)
set.seed(20260227)

ph_true <- ph(
  structure = "merlang",
  block_sizes = c(4,4),
  probs = c(0.25, 0.75),
  rates = c(0.1, 0.02)
)

n <- 1000
y <- sim(ph_true, n = n)

hist(y,probability = T,breaks=100)
sq <- seq(0,1000,0.1)
lines(sq,dens(ph_true,sq))

ft <- fit(
  ph(structure = "merlang", block_sizes = c(4,4)),
  y,
  stepsEM = 400
)

lines(sq,dens(ft,sq),col="red")

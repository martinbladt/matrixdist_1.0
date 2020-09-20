# MPH
alpha <- c(0.15, 0.85)
pi <- c(alpha, 0 ,0)
T11 <- matrix(c(c(-2,9),c(0,-11)), nrow = 2, ncol = 2)
T12 <- matrix(c(c(2,0),c(0,2)), nrow = 2, ncol = 2)
T22 <- matrix(c(c(-1,0),c(0.5,-5)), nrow = 2, ncol = 2)
T <- merge_matrices(T11, T12, T22)
R <- matrix(c(c(1,1,0,0), c(0,0,1,1)), ncol=2)
n <- 100
set.seed(1)
x <- rmph(n, pi, T, R) 


xweight <- rep(1, length(x[,1]))


mph_par <- random_structure(6)
R_rnd <- random_reward(6,2)

pi_fit <- clone_vector(mph_par$pi)
T_fit <- clone_matrix(mph_par$T)
R_fit <-clone_matrix(R_rnd)

x_sum <- sum_data(x)

logLikelihood <- function(x, pi, T) {
  sum(log(phdensity(x,pi,T)));
}

logLikelihood(x, pi_fit, T_fit)

censored <- numeric(0)
rcweight <- numeric(0)


stepsEM1 <- 100

for (k in 1:stepsEM) {
  
  EMstep(pi_fit, T_fit, x_sum, xweight, censored, rcweight)
  
  if (k %% 10 == 0) {
    print(logLikelihood(x_sum, pi_fit, T_fit))
  }
}
pi_fit
T_fit
R_fit

library(progress)

stepsEM2 <- 100

pb <- progress_bar$new(total = stepsEM2)

for (k in 1:stepsEM) {
  
  secondEMstep(x, xweight, matrix(censored), rcweight, pi_fit, T_fit, R_fit)
  
  pb$tick()
}
pi_fit
T_fit
R_fit




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


logLikelihood <- function(x, alpha, T11, T12, T22) {
  sum(log(bivphden(x,alpha,T11,T12,T22)));
}


logLikelihood(x, alpha, T11, T12, T22)


mph_par <- random_phase_BivPH(3,3)

alpha_fit <- clone_vector(mph_par$alpha)
T11_fit <- clone_matrix(mph_par$T11)
T12_fit <- clone_matrix(mph_par$T12)
T22_fit <- clone_matrix(mph_par$T22)


logLikelihood(x, alpha_fit, T11_fit, T12_fit, T22_fit)

stepsEM <- 1000

for (k in 1:stepsEM) {
  
  EMstep_bivph(x, xweight, alpha_fit, T11_fit, T12_fit, T22_fit)
  
  if (k %% 10 == 0) {
    print(logLikelihood(x, alpha_fit, T11_fit, T12_fit, T22_fit))
  }
}
alpha_fit
T11_fit
T12_fit
T22_fit

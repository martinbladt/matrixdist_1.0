# Bivariate PH
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


object <- bph(alpha, T11, T12, T22)
dens(object, x)

data <- sim(object)

F1 <- fit(object, data, stepsEM = 200)

obj2 <- ibph(object, gfun = "Pareto", gfun_pars = c(2, 1))
data <- sim(obj2)

F2 <- fit(obj2, data, stepsEM = 100)

x <- obj2
y <- data

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



#Bivariate inhomogenous PH
alpha <- c(0.15, 0.85)
pi <- c(alpha, 0 ,0)
T11 <- matrix(c(c(-2,9),c(0,-11)), nrow = 2, ncol = 2)
T12 <- matrix(c(c(2,0),c(0,2)), nrow = 2, ncol = 2)
T22 <- matrix(c(c(-1,0),c(0.5,-5)), nrow = 2, ncol = 2)
T <- merge_matrices(T11, T12, T22)
R <- matrix(c(c(1,1,0,0), c(0,0,1,1)), ncol=2)
beta <- c(0.4, 0.7)
n <- 100
set.seed(1)
x <- rimph(n, "Weibull", pi, T, R, beta) 

object <- bph(alpha, T11, T12, T22)
obj2 <- ibph(object, gfun = "Weibull", gfun_pars = beta)
data <- x
F2 <- fit(obj2, data, stepsEM = 100)

xweight <- rep(1, length(x[,1]))


logLikelihoodIPH <- function(x, alpha, T11, T12, T22, beta) {
  return(sum(log(bivmWeibullden(x,alpha,T11,T12,T22, beta))))
}


mLL <- function(x, alpha, T11, T12, T22, beta) {
  return(-logLikelihoodIPH(x, alpha, T11, T12, T22, beta))
}


logLikelihoodIPH(x, alpha, T11, T12, T22, beta)


mph_par <- random_phase_BivPH(3,3)

alpha_fit <- clone_vector(mph_par$alpha)
T11_fit <- clone_matrix(mph_par$T11)
T12_fit <- clone_matrix(mph_par$T12)
T22_fit <- clone_matrix(mph_par$T22)
beta_fit <- c(0.2, 0.6)

logLikelihoodIPH(x, alpha_fit, T11_fit, T12_fit, T22_fit, beta_fit)

x_transform <- x

stepsEM <- 100

for (k in 1:stepsEM) {
  x_transform[,1] <- x[,1]^(beta_fit[1])
  x_transform[,2] <- x[,2]^(beta_fit[2])
  
  EMstep_bivph(x_transform, xweight, alpha_fit, T11_fit, T12_fit, T22_fit)
  
  opt <- suppressWarnings(optim(par = beta_fit, fn = mLL,  x = x, alpha = alpha_fit, T11 = T11_fit, T12 = T12_fit, T22 = T22_fit))
  beta_fit <- opt$par
  
  if (k %% 10 == 0) {
    print(- opt$value)
  }
}
alpha_fit
T11_fit
T12_fit
T22_fit
beta_fit

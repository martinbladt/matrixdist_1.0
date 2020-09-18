library(bench)

alpha <- c(0.5, 0.3, 0.2)
T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
beta <- 0.5 
n <- 1000
set.seed(1)
x <- riph(n, "Weibull", alpha, T, beta) 

xweight <- rep(1, length(x))

observations <- cbind(x, xweight)
mat <- data.frame(observations)
names(mat) <- c('obs', 'weight')

# aggregate(mat$weight, by = list(obs = mat$obs), FUN=sum)

sort_obs <- sort(unique(x))

cum_weight <- NULL
for(i in sort_obs){
  cum_weight <- c(cum_weight, sum(mat$weight[which(mat$obs==i)]))
}

right_cen <- numeric(0)
right_cen_weight <- numeric(0)


loglikelihoodMW <- function(beta, alpha, T, x) {
  return(sum(log(mweibullden(x, alpha, T, beta))))
}

RKstep <- default_step_length(T)


loglikelihoodMW(beta, alpha, T, x)

g <- function(x, beta) { x^(1/beta) }
g_inv <- function(x, beta) { x^beta}
lambda <- function(x, beta) {beta * x^(beta - 1)}

loglikelihoodMW2 <- function(beta, alpha, T, x) {
  return(sum(log(iphdensity(x, alpha, T, g, g_inv, lambda, beta))))
}

logLikelihoodIPH_RK(RKstep, alpha, T, g, g_inv, lambda, beta, sort_obs,cum_weight, right_cen, right_cen_weight) 

logLikelihoodMWeib_RK(RKstep, alpha, T, beta, sort_obs, cum_weight, right_cen, right_cen_weight) 


bench::mark(
  loglikelihoodMW(beta, alpha, T, x),
  loglikelihoodMW2(beta, alpha, T, x),
  logLikelihoodIPH_RK(RKstep, alpha, T, g, g_inv, lambda, beta, sort_obs,cum_weight, right_cen, right_cen_weight),
  logLikelihoodMWeib_RK(RKstep, alpha, T, beta, sort_obs, cum_weight, right_cen, right_cen_weight) 
)[1:6]


logLikelihoodMPar_RK(RKstep, alpha, T, beta, sort_obs, cum_weight, right_cen, right_cen_weight) 
sum(log(mparetoden(x, alpha, T, beta)))

betavec <- c(1, 2, 0)
logLikelihoodMGEV_RK(RKstep, alpha, T, betavec, sort_obs, cum_weight, right_cen, right_cen_weight) 
sum(log(mGEVden(x, alpha, T, 1, 2, 0)))


x2 <- riph(n, "Gompertz", alpha, T, beta) 
sort_obs2 <- sort(unique(x2))
logLikelihoodMGomp_RK(RKstep, alpha, T, beta, sort_obs2, cum_weight, right_cen, right_cen_weight) 
sum(log(mgompertzden(x2, alpha, T, beta)))



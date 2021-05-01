# Mixing component
density <- function(t, parameters) {
  return(dgamma(t, shape = parameters[1], rate = parameters[2]))
}

parameters <- c(1,2)
x <- 1
density(x, parameters)


# Discretization of the continuous distribution for calculations
truncationPoint = 10
maxProbability = 0.01
maxDeltat = 0.5

gridDensity <- function(density, parameters, truncationPoint, maxProbability, maxDeltat) {
  deltat <- 0
  t <- 0.0001 # It should work with 0
  
  prob = numeric(0)
  value = numeric(0)
  
  j <- 1
  
  while (t < truncationPoint) {
    if (density(t, parameters) < maxProbability / maxDeltat) {
      deltat = maxDeltat
    }
    else {
      deltat = maxProbability / density(t, parameters)
    }
    proba_aux = deltat / 6 * (density(t, parameters) + 4 * density(t + deltat / 2, parameters) + density(t + deltat, parameters))
    while (proba_aux > maxProbability) {
      deltat = deltat * 0.9
      proba_aux = deltat / 6 * (density(t ,parameters) + 4 * density(t + deltat / 2, parameters) + density(t + deltat, parameters))
    }
    if (proba_aux > 0) {
      value[j] = (t * density(t, parameters) + 4 * (t + deltat / 2) * density(t + deltat / 2, parameters) + (t + deltat) * density(t + deltat, parameters)) / (density(t, parameters) + 4 * density(t + deltat / 2, parameters) + density(t + deltat, parameters))
      prob[j] = proba_aux
      j <- j + 1
    }
    t = t + deltat
  }
  prob <- prob/sum(prob)
  my_list <- list("val" = value, "prob" = prob)
  return(my_list) 
  
} # se deberia poder hacer esto con la cdf?


grid <- gridDensity(density, parameters, truncationPoint, maxProbability, maxDeltat) 


base::plot(grid$val, grid$prob)

sum(grid$prob)

#Initial values of the distribution be fitted
ph1 <- ph(structure = "GCoxian", dimension = 4)
pi <- ph1@pars$alpha
S <- ph1@pars$S
I <- diag(4)
alpha <- 1
beta <- 2 
parameters <- c(1,2)

# Density of the mPH
library(expm)
computeDen <- function(theSample, pi, S, parameters) {
  den = numeric(0)
  for (i in 1:length(theSample)) {
    den[i] = parameters[1] * parameters[2]^parameters[1] * pi %*% expm( -(parameters[1] + 1) * logm(parameters[2] * I - theSample[i] * S ) ) %*% (S * (-1)) %*% rep(1,length(pi))
  }
  return(den)
}

computeDen2 <- function(theSample, pi, S, parameters) {
  den = numeric(0)
  for (i in 1:length(theSample)) {
    eg <- eigen(parameters[2] * I - theSample[i] * S)
    p <- eg$vectors
    d <- diag(eg$values^{-(parameters[1] + 1)})
    mm <- p %*% d %*% solve(p)
    den[i] = parameters[1] * parameters[2]^parameters[1] * pi %*%  mm %*% (S * (-1)) %*% rep(1,length(pi))
  }
  return(den)
}




computeTail <- function(theSample, pi, S, parameters) {
  den = numeric(0)
  for (i in 1:length(theSample)) {
    den[i] = parameters[2]^parameters[1] * pi %*% expm( -(parameters[1]) * logm(parameters[2] * I - theSample[i] * S ) ) %*% rep(1,length(pi))
  }
  return(den)
}

thePar <- c(alpha, beta)
computeDen(theSample, pi, S, thePar) 
computeDen2(theSample, pi, S, thePar) 


# Sample of a Pareto as a test
library(actuar)
theSample <- rpareto(1000, shape = 2.5, scale = 1)
sum(log(dpareto(theSample,shape = 2.5, scale = 1)))
n <- length(theSample)

xweight <- rep(1, length(theSample))
observations <- cbind(theSample, xweight)
mat <- data.frame(observations)
names(mat) <- c('obs', 'weight')

sort_obs <- sort(unique(theSample))

cum_weight <- NULL
for(i in sort_obs){
  cum_weight <- c(cum_weight, sum(mat$weight[which(mat$obs==i)]))
} #need to be sorted still... use data aggregation function


# Function needed in the M step for theta
maxfunction <- function(parameters,  density, valfn, L) {
  cum <- 0
  for (i in 1:(length(L))) {
    cum = cum - log(density(valfn[i], parameters)) * L[i]
  }
  return(cum)
}


library(matrixdist)
# EM 
stepsEM <- 500
t1 <- t2 <- t3 <- t4 <- t5 <- 0
for (k in 1:stepsEM) {
  RKstep <- default_step_length(S)

  theDensity <- computeDen2(sort_obs, pi, S, parameters) 

  grid <- gridDensity(density, parameters, truncationPoint, maxProbability, maxDeltat) 

  L <- EMstep_mPH_RK(RKstep, pi, S, sort_obs, cum_weight,  theDensity, grid$val, grid$prob)

  pi <- pi/sum(pi) #To compensate for the error due  to discretization - double check
  
  parameters <- suppressWarnings(optim(par = parameters, fn = maxfunction, density= density, valfn = grid$val, L = L, list(maxit = 10))$par)

  if (k %% 10 == 0) {
    print(c(k, sum(log(theDensity))))
  }
  
}

theDensity <- computeDen2(sort_obs, pi, S, parameters) 
hist(theSample, breaks = 40 ,freq = FALSE)
lines(sort_obs, theDensity, col = "blue")

base::plot(sort_obs, log(theDensity), type = "l", col = "blue")
lines(sort_obs, log(dpareto(sort_obs,shape = 2.5, scale = 1)), col = "red")

#likelihood better than using original distribution


## ToDo 
# I saw the likelihood slighly decrease in some cases / check why / due to discretization?
# See if I can improve speed by using the derivative / this will depend on the mixing distribution
# qq-plot would be better to see quality of fit
# Simulation of this type of RV, it should be easy
# Compare with a matrix-Pareto type I?

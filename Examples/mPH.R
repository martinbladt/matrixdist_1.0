density <- function(t, parameters) {
  return(dgamma(t, shape = parameters[1], rate = parameters[2]))
}

parameters <- c(1,2)
x <- 1
density(x, parameters)

truncationPoint = 10

maxProbability = 0.1

maxDeltat = 0.1


gridDensity <- function(density, parameters, truncationPoint, maxProbability, maxDeltat) {
  deltat <- 0
  t <- 0
  
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
  
  my_list <- list("val" = value, "prob" = prob)
  return(my_list) 
  
}


grid <- gridDensity(density, parameters, truncationPoint, maxProbability, maxDeltat) 

library(actuar)
theSample <- rpareto(100, shape = 2.5, scale = 1)
n <- length(theSample)

library(expm)
pi <- c(0.15, 0.85)
S <- matrix(c(c(-1,0.5),c(0,-5)), nrow = 2, ncol = 2)
I <- matrix(c(c(1,0),c(0,1)), nrow = 2, ncol = 2)
alpha <- 1
beta <- 2 

x <- 1
alpha * beta^alpha * pi %*% expm( -(alpha + 1) * logm(beta * I - x * S ) ) %*% (S * (-1)) %*% rep(1,length(pi))


computeDen <- function(theSample, pi, S, parameters) {
  den = numeric(0)
  for (i in 1:length(theSample)) {
    den[i] = parameters[1] * parameters[2]^parameters[1] * pi %*% expm( -(parameters[1] + 1) * logm(parameters[2] * I - theSample[i] * S ) ) %*% (S * (-1)) %*% rep(1,length(pi))
  }
  return(den)
}

thePar <- c(alpha, beta)
computeDen(theSample, pi, S, thePar) 




stepsEM <- 1000

for (k in 1:stepsEM) {
  RKstep <- default_step_length(S)
  
  L <- EMstep_mPH_RK(RKstep, pi, S, obs, weight,  density, valmix, probmix)
  
  
}



ph1 <- ph(structure = "GCoxian", dimension = 5)
ph1 <- ph(structure = "Coxian", dimension = 4)
set.seed(1)
ph1 <- ph(structure = "General", dimension = 4)
coef(ph1)
ph_par <- ph1@pars
alpha <- clone_vector(ph_par$alpha)
S <- clone_matrix(ph_par$S)
theta <- 2

x <- 1
mp3density(x, 1, alpha, S)

phLaplace(x, alpha, S) 

phdensity(x, alpha, S)


mp3cdf(10, 1, alpha, S)

laplaceph3(1, alpha, S)


loglikelihood_mp3 <- function(y, theta, alpha, S) {
  sum(log(mp3density(y, theta, alpha, S)))
}

loglikelihood_mp3_cens <- function(y, cens, theta, alpha, S) {
  sum(log(mp3density(y, theta, alpha, S))) + sum(log(1 - mp3cdf(cens, theta, alpha, S)))
}


Ezgiveny <- function(thetamax, y, theta, alpha, S) {
  -sum(log(thetamax * y^(thetamax-1))  - y^thetamax * 2 * theta * y^(theta-1) * laplaceph3(y^theta, alpha, S) / mp3density(y, theta, alpha, S))
}


Ezgivenycen <- function(thetamax, y, cens, theta, alpha, S) {
  -sum(log(thetamax * y^(thetamax-1))  - y^thetamax * 2 * theta * y^(theta-1) * laplaceph3(y^theta, alpha, S) / mp3density(y, theta, alpha, S)) + sum(cens^thetamax * laplaceph2(cens^theta, alpha, S) / (1- mp3cdf(cens, theta, alpha, S)) )
}

mp3density(y, theta, alpha, S)

y <- rgamma(500, 2, 1)
cens <- rgamma(500, 2, 1)
Ezgiveny(theta, y,  theta, alpha, S) 
Ezgivenycen(theta, y, cens, theta, alpha, S) 

theta_fit <- optim(par = theta, fn = Ezgiveny, y = y, theta = theta, alpha = alpha, S = S)$par

fn(theta)
optim(par = theta, fn)$par


conditional_density <- function(z, y, theta, alpha, S) {
  mean(z * theta * y^(theta-1) * exp(- z * y^theta) * phdensity(z, alpha, S) / mp3density(y, theta, alpha, S))
}
conditional_density(10, y, theta, alpha, S)


conditional_density_cens <- function(z, y, cens, theta, alpha, S) {
  (sum(z * theta * y^(theta-1) * exp(- z * y^theta) * phdensity(z, alpha, S) / mp3density(y, theta, alpha, S)) + sum(exp(- z * cens^theta) * phdensity(z, alpha, S) / (1 - mp3cdf(cens, theta, alpha, S)))) / (length(y) + length(cens))
}
conditional_density_cens(10, y,cens, theta, alpha, S)

condendity(10, y, alpha, S, theta)




#Discretization of density
truncationPoint <- 5
maxProbability <-  0.01
maxDeltat <- 0.05

deltat <- 0
t <- 0.0001 

prob = numeric(0)
value = numeric(0)

j <- 1

while (t < truncationPoint) {
  if (conditional_density(t, y, theta, alpha, S) < maxProbability / maxDeltat) {
    deltat = maxDeltat
  }
  else {
    deltat = maxProbability / conditional_density(t, y, theta, alpha, S)
  }
  proba_aux = deltat / 6 * (conditional_density(t, y, theta, alpha, S) + 4 * conditional_density(t + deltat / 2, y, theta, alpha, S) + conditional_density(t + deltat, y, theta, alpha, S) )
  while (proba_aux > maxProbability) {
    deltat = deltat * 0.9
    proba_aux = deltat / 6 * (conditional_density(t, y, theta, alpha, S) + 4 * conditional_density(t + deltat / 2, y, theta, alpha, S) + conditional_density(t + deltat, y, theta, alpha, S))
  }
  if (proba_aux > 0) {
    value[j] = (t * conditional_density(t, y, theta, alpha, S)  + 4 * (t + deltat / 2) * conditional_density(t + deltat / 2, y, theta, alpha, S) + (t + deltat) * conditional_density(t + deltat, y, theta, alpha, S)) / (conditional_density(t, y, theta, alpha, S) + 4 * conditional_density(t + deltat / 2, y, theta, alpha, S) + conditional_density(t + deltat, y, theta, alpha, S))
    prob[j] = proba_aux
    j <- j + 1
  }
  t = t + deltat
}

length(prob)
sum(prob)


data <- read.table("/Users/jorgeyslas/Documents/Research/Frailty/danish.txt", header = TRUE)
claims <- data$total
evir::hill(y)
y <- claims[claims>1] - 1


loglikelihood_mp3(y, theta, alpha, S)

dt = 0.01
thez <- seq(0.01, 10, by = dt)
weights <- rep(0, length(thez))

NumSteps <- 40

for (k in 1:NumSteps) {
  
  theta_fit <- suppressWarnings(optim(par = theta, fn = Ezgiveny, y = y, theta = theta, alpha = alpha, S = S)$par)
  
  # for (l in 1:length(thez)) {
  #  weights[l] <- dt * conditional_density(thez[l], y, theta, alpha, S)
  # }
  #Discretization of density
  truncationPoint <- 20
  maxProbability <-  0.0025
  maxDeltat <- 0.05

  deltat <- 0
  t <- 0.0001

  prob = numeric(0)
  value = numeric(0)

  j <- 1

  while (t < truncationPoint) {
    if (conditional_density(t, y, theta, alpha, S) < maxProbability / maxDeltat) {
      deltat = maxDeltat
    }
    else {
      deltat = maxProbability / conditional_density(t, y, theta, alpha, S)
    }
    proba_aux = deltat / 6 * (conditional_density(t, y, theta, alpha, S) + 4 * conditional_density(t + deltat / 2, y, theta, alpha, S) + conditional_density(t + deltat, y, theta, alpha, S) )
    while (proba_aux > maxProbability) {
      deltat = deltat * 0.9
      proba_aux = deltat / 6 * (conditional_density(t, y, theta, alpha, S) + 4 * conditional_density(t + deltat / 2, y, theta, alpha, S) + conditional_density(t + deltat, y, theta, alpha, S))
    }
    if (proba_aux > 0) {
      value[j] = (t * conditional_density(t, y, theta, alpha, S)  + 4 * (t + deltat / 2) * conditional_density(t + deltat / 2, y, theta, alpha, S) + (t + deltat) * conditional_density(t + deltat, y, theta, alpha, S)) / (conditional_density(t, y, theta, alpha, S) + 4 * conditional_density(t + deltat / 2, y, theta, alpha, S) + conditional_density(t + deltat, y, theta, alpha, S))
      prob[j] = proba_aux
      j <- j + 1
    }
    t = t + deltat
  }
 
  auxph <- ph(alpha = alpha, S = S)
  auxfit <- fit(x = auxph, y = value, weight = prob, stepsEM = 99)
  #auxfit <- fit(x = auxph, y = thez, weight = weights, stepsEM = 50)
  
  alpha <- clone_vector(coef(auxfit)$alpha)
  S <- clone_matrix(coef(auxfit)$S) 
  theta <- theta_fit
  
  cat("\r", "iteration:", k, ", logLik:", loglikelihood_mp3(y, theta, alpha, S), ", theta:" , theta, ", sumweights:" , sum(prob), sep = " ")
  
}
S

sq <- seq(0.01, 10, by = 0.01)
hist(y, breaks = 1800, freq = FALSE, xlim = c(0,8), xlab = "x" , ylab = NULL, main = "Histogram vs fitted density",  cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
lines(sq, mp3density(sq, theta, alpha, S), col = "#3498DB", lwd = 2.5)
legend("topright", 
       legend = "Matrix-Pareto type III ", 
       col ="#3498DB", 
       bty = "n", 
       lwd = 2.5, 
       cex = 1.1, 
       text.col = "black", 
       horiz = FALSE , 
       inset = c(0.05))



# QQ-plot
s <- 50000
qarray <- array(1:s)
x0 <- 0
delta <- max(y) * 1.5 / s
for(k in 1:s) {
  qarray[k] <- mp3cdf(x0 + k * delta, theta, alpha, S)
}
quant <- array(1:199)

for (k in 1:199) {
  quant[k] <- x0 + (min(which(qarray > k / 200.0))) * delta
}

samplequant <- quantile(y, probs = seq(0.005, 1 - 0.005, 0.005))


lambda <- 1/mean(log(y +1))
quant <- (1-seq(0.005, 1 - 0.005, 0.005))^(-1/lambda) - 1

base::plot( samplequant, quant,xlab = "Sample", ylab = "Fitted distribution", main = "QQ-plot - Matrix-Pareto type III",  cex.main=1.8, cex.lab=1.5, cex.axis=1.5, bty="l")
abline(0,1)




############################
# Loss data
############################

# scaled using 10,000 as in the CPH paper
scale_factor <- 10000
data(loss, package="copula")
y <- loss$loss[loss$censored == 0] / scale_factor
rcens <- loss$loss[loss$censored == 1] / scale_factor


# ph1 <- ph(structure = "Coxian", dimension = 3)
# NumSteps <- 30
# With simple discretization

set.seed(1)
ph1 <- ph(structure = "Coxian", dimension = 4)
coef(ph1)
ph_par <- ph1@pars
alpha <- clone_vector(ph_par$alpha)
S <- clone_matrix(ph_par$S)
theta <- 2



loglikelihood_mp3_cens(y, rcens, theta, alpha, S)

dt = 0.01
thez <- seq(0.01, 15, by = dt)
weights <- rep(0, length(thez))

NumSteps <- 40

for (k in 1:NumSteps) {
  
  theta_fit <- suppressWarnings(optim(par = theta, fn = Ezgivenycen, y = y, cens = rcens, theta = theta, alpha = alpha, S = S)$par)
  
  # for (l in 1:length(thez)) {
  #   weights[l] <- dt * conditional_density_cens(thez[l], y, rcens, theta, alpha, S)
  # }
  #Discretization of density
  truncationPoint <- 18
  maxProbability <-  0.005
  maxDeltat <- 0.05

  deltat <- 0
  t <- 0.0001

  prob = numeric(0)
  value = numeric(0)

  j <- 1

  while (t < truncationPoint) {
    if (conditional_density_cens(t, y, rcens, theta, alpha, S) < maxProbability / maxDeltat) {
      deltat = maxDeltat
    }
    else {
      deltat = maxProbability / conditional_density_cens(t, y, rcens, theta, alpha, S)
    }
    proba_aux = deltat / 6 * (conditional_density_cens(t, y, rcens, theta, alpha, S) + 4 * conditional_density_cens(t + deltat / 2, y, rcens, theta, alpha, S) + conditional_density_cens(t + deltat, y, rcens, theta, alpha, S) )
    while (proba_aux > maxProbability) {
      deltat = deltat * 0.9
      proba_aux = deltat / 6 * (conditional_density_cens(t, y, rcens, theta, alpha, S) + 4 * conditional_density_cens(t + deltat / 2, y, rcens, theta, alpha, S) + conditional_density_cens(t + deltat, y, rcens, theta, alpha, S))
    }
    if (proba_aux > 0) {
      value[j] = (t * conditional_density_cens(t, y, rcens, theta, alpha, S)  + 4 * (t + deltat / 2) * conditional_density_cens(t + deltat / 2, y, rcens, theta, alpha, S) + (t + deltat) * conditional_density_cens(t + deltat, y, rcens, theta, alpha, S)) / (conditional_density_cens(t, y, rcens, theta, alpha, S) + 4 * conditional_density_cens(t + deltat / 2, y, rcens, theta, alpha, S) + conditional_density_cens(t + deltat, y, rcens, theta, alpha, S))
      prob[j] = proba_aux
      j <- j + 1
    }
    t = t + deltat
  }

  auxph <- ph(alpha = alpha, S = S)
  auxfit <- fit(x = auxph, y = value, weight = prob, stepsEM = 99)
  #auxfit <- fit(x = auxph, y = thez, weight = weights, stepsEM = 99)

  alpha <- clone_vector(coef(auxfit)$alpha)
  S <- clone_matrix(coef(auxfit)$S)
  theta <- theta_fit

  cat("\r", "iteration:", k, ", logLik:", loglikelihood_mp3_cens(y, rcens, theta, alpha, S), ", theta:" , theta, ", sumweights:" , sum(prob), sep = " ")

}

round(S, digits = 4)
round(theta, digits = 4)


library(survival)
lossmod <- loss
lossmod$loss <- lossmod$loss / scale_factor
sq <- seq(0, 2500, by = 0.1)
surv <- survfit(Surv(loss, 1 - censored) ~ 1, data = lossmod)
png("loss_fit.png", width = 8, height = 8, units = 'in', res = 200)
base::plot(surv$time ,surv$cumhaz, 
           xlab = "Loss", 
           ylab = "Cumulative hazard",
           main = "Fit to loss data", type = "l",  cex.main=1.8, cex.lab=1.5, cex.axis=1.5, bty="l")
points(surv$time ,surv$cumhaz, pch = 4, cex = 0.5)
lines(sq, -log( 1- mp3cdf(sq, theta, alpha, S)), col = "#3498DB", lwd = 2.5)
legend("bottomright", 
       legend = "Matrix-Pareto type III", 
       col ="#3498DB", 
       bty = "n", 
       lwd = 2.5, 
       cex = 1.1, 
       text.col = "black", 
       horiz = FALSE , 
       inset = c(0.05))
dev.off()



############################
# lognormal fraitly
############################
n1 <- 500

set.seed(1)
z <- rlnorm(n1, -0.35, 0.8)

# Gompertz 
b = 0.01 
c = 1

group1 <- (1/c) * log(1 - c * log(runif(n1))/(z * b))

n2 <- 500
beta <- 0.5
set.seed(12)
z2 <- rlnorm(n2, -0.35, 0.8) * exp(beta)

group2 <- (1/c) * log(1 - c * log(runif(n1))/(z2 * b))


y0 <- c(group1, group2) 
x0 <- c(rep(0, n1), rep(1, n2))

hist(y0)

# here the par_max is of the form (b,c,beta)
Ezgivenycov <- function(par_max, y, x, theta, alpha, S) {
  -sum( x * par_max[3] + log(par_max[1])  + y * par_max[2]  - exp(x * par_max[3]) * par_max[1] * (exp(par_max[2] * y) - 1) / par_max[2] * 2 * laplaceph3(exp(x * theta[3]) * theta[1] * (exp(theta[2] * y) -1 ) / theta[2], alpha, S) / laplaceph2(exp(x * theta[3]) * theta[1] * (exp(theta[2] * y) -1 ) / theta[2], alpha, S))
}

set.seed(1)
phl <- ph(structure = "General", dimension = 2)
coef(phl)
ph_par <- phl@pars
alpha <- clone_vector(ph_par$alpha)
S <- clone_matrix(ph_par$S)
theta <- c(0.01,1,0.5)

base::plot(y0, laplaceph3(exp(x0 * theta[3]) * theta[1] * (exp(theta[2] * y0) -1 ) / theta[2], alpha, S) / laplaceph2(exp(x0 * theta[3]) * theta[1] * (exp(theta[2] * y0) -1 ) / theta[2], alpha, S))

lnphden <- function(y, x, theta, alpha, S) {
  exp(theta[3] * x) * theta[1] * exp(theta[2] * y ) * laplaceph2(exp(x * theta[3]) * theta[1] * (exp(theta[2] * y) -1 ) / theta[2], alpha, S)
}

sq <- seq(0, 10, by = 0.001)
integral <- 0
for(i in 1:length(sq)){
  integral <- integral + lnphden(sq[i], 1, theta, alpha, S)
}
0.001 * integral 

lnphden(1, 1, theta, alpha, S)

loglikelihood_ln <- function(y, x, theta, alpha, S) {
  sum( log(exp(theta[3] * x) * theta[1] * exp(theta[2] * y ) * laplaceph2(exp(x * theta[3]) * theta[1] * (exp(theta[2] * y) -1 ) / theta[2], alpha, S)))
}

library(lamW)
lambertW0(exp(1))

laplacenormal0 <- function(s, sigma) {
  (1 / sqrt(1 + lambertW0(s * sigma^2))) * exp(-(0.5/sigma^2 ) * (lambertW0(s * sigma^2))^2 - (1 /sigma^2)  * lambertW0(s * sigma^2) )
}

laplacenormal0(0.5,1)


laplacenormal <- function(s, mu, sigma) {
  laplacenormal0(s * exp( mu), sigma) 
}

sum(log(exp(x0 * beta) * b * exp(c * y0) * laplacenormal( exp(x0 * beta) * b * (exp(c * y0 ) - 1) / c,  -0.35, 0.8)))


loglikelihood_ln(y0, x0, theta, alpha, S) 

Ezgivenycov(c(1,0.01,1), y0, x0, theta, alpha, S)

par_ini <- c(0.01,1,0.5)

par_fit <- optim(par = par_ini, fn = Ezgivenycov, y = y0, x = x0 , theta = theta, alpha = alpha, S = S)$par

loglikelihood_ln(y0, x0, par_fit, alpha, S) 



conditional_density_ln <- function(z, y, x, theta, alpha, S) {
  mean( z * exp(theta[3] * x) * theta[1] * exp(theta[2] * y) * exp(- z * exp(theta[3] * x) * theta[1] * (exp(theta[2] * y) - 1) / theta[2] ) * phdensity(z, alpha, S) / lnphden(y, x, theta, alpha, S) )
}
conditional_density_ln(10, y0, x0, theta, alpha, S)
lnphden(10, 1, theta, alpha, S)

sq <- seq(0, 10, by = 0.01)
integral <- 0
for(i in 1:length(sq)){
  integral <- integral + conditional_density_ln(sq[i], y0, x0, theta, alpha, S)
}
0.01 * integral 



loglikelihood_ln(y0, x0, theta, alpha, S) 

dt = 0.01
thez <- seq(0.01, 15, by = dt)
weights <- rep(0, length(thez))

NumSteps <- 100

for (k in 1:NumSteps) {
  
  theta_fit <- suppressWarnings(optim(par = theta, fn = Ezgivenycov, y = y0, x = x0, theta = theta, alpha = alpha, S = S)$par)
  
  # for (l in 1:length(thez)) {
  #   weights[l] <- dt * conditional_density_cens(thez[l], y, rcens, theta, alpha, S)
  # }
  #Discretization of density
  truncationPoint <- 10
  maxProbability <-  0.005
  maxDeltat <- 0.05
  
  deltat <- 0
  t <- 0.0001
  
  prob = numeric(0)
  value = numeric(0)
  
  j <- 1
  
  while (t < truncationPoint) {
    if (conditional_density_ln(t, y0, x0, theta, alpha, S) < maxProbability / maxDeltat) {
      deltat = maxDeltat
    }
    else {
      deltat = maxProbability / conditional_density_ln(t, y0, x0, theta, alpha, S)
    }
    proba_aux = deltat / 6 * (conditional_density_ln(t, y0, x0, theta, alpha, S) + 4 * conditional_density_ln(t + deltat / 2, y0, x0, theta, alpha, S) + conditional_density_ln(t + deltat, y0, x0, theta, alpha, S) )
    while (proba_aux > maxProbability) {
      deltat = deltat * 0.9
      proba_aux = deltat / 6 * (conditional_density_ln(t, y0, x0, theta, alpha, S) + 4 * conditional_density_ln(t + deltat / 2, y0, x0, theta, alpha, S) + conditional_density_ln(t + deltat, y0, x0, theta, alpha, S))
    }
    if (proba_aux > 0) {
      value[j] = (t * conditional_density_ln(t, y0, x0, theta, alpha, S)  + 4 * (t + deltat / 2) * conditional_density_ln(t + deltat / 2, y0, x0, theta, alpha, S) + (t + deltat) * conditional_density_ln(t + deltat, y0, x0, theta, alpha, S)) / (conditional_density_ln(t, y0, x0, theta, alpha, S) + 4 * conditional_density_ln(t + deltat / 2, y0, x0, theta, alpha, S) + conditional_density_ln(t + deltat, y0, x0, theta, alpha, S))
      prob[j] = proba_aux
      j <- j + 1
    }
    t = t + deltat
  }
  
  auxph <- ph(alpha = alpha, S = S)
  auxfit <- fit(x = auxph, y = value, weight = prob, stepsEM = 99)
  #auxfit <- fit(x = auxph, y = thez, weight = weights, stepsEM = 99)
  
  alpha <- clone_vector(coef(auxfit)$alpha)
  S <- clone_matrix(coef(auxfit)$S)
  theta <- theta_fit
  
  cat("\r", "iteration:", k, ", logLik:", loglikelihood_ln(y0, x0, theta, alpha, S), ", theta:" , theta, ", sumweights:" , sum(prob), sep = " ")
  
}

round(alpha, digits = 4)
round(S, digits = 4)
round(theta, digits = 4)


lnphcdf <- function(y, x, theta, alpha, S) {
  1 - phLaplace(exp(x * theta[3]) * theta[1] * (exp(theta[2] * y) -1 ) / theta[2], alpha, S)
}
lnphcdf(10, 0, theta, alpha, S)


# QQ-plot
s <- 50000
qarray0 <- array(1:s)
qarray1 <- array(1:s)
x0 <- 0
delta <- max(y0) * 1.5 / s
for(k in 1:s) {
  qarray0[k] <- lnphcdf(x0 + k * delta, 0, theta, alpha, S)
  qarray1[k] <- lnphcdf(x0 + k * delta, 1, theta, alpha, S)
}
quant0 <- array(1:199)
quant1 <- array(1:199)

for (k in 1:199) {
  quant0[k] <- x0 + (min(which(qarray0 > k / 200.0))) * delta
  quant1[k] <- x0 + (min(which(qarray1 > k / 200.0))) * delta
}

samplequant0 <- quantile(group1, probs = seq(0.005, 1 - 0.005, 0.005))
samplequant1 <- quantile(group2, probs = seq(0.005, 1 - 0.005, 0.005))

library(latex2exp)
png("ln_qq0.png", width = 8, height = 8, units = 'in', res = 200)
base::plot( samplequant0, quant0,xlab = TeX("Sample - $X = 0$ "), ylab = "Fitted distribution", main = "QQ-plot - PH frailty",  cex.main=1.8, cex.lab=1.5, cex.axis=1.5, bty="l")
abline(0,1)
dev.off()

png("ln_qq1.png", width = 8, height = 8, units = 'in', res = 200)
base::plot( samplequant1, quant1,xlab = TeX("Sample - $X = 1$ "), ylab = "Fitted distribution", main = "QQ-plot - PH frailty",  cex.main=1.8, cex.lab=1.5, cex.axis=1.5, bty="l")
abline(0,1)
dev.off()

# Histograms
hist(group1, freq = F, breaks = 30, xlab = NULL, main = TeX("\\textbf{Histogram -} $\\beta = 0$ "))
q <- seq(0.001, 9, by = 0.01)
lines(q, lnphden(q, rep(0,length(q) ), theta, alpha, S), col = "#3498DB", lwd = 2.5)

hist(group2, freq = F, breaks = 30, xlab = NULL, main = TeX("\\textbf{Histogram -} $\\beta = 1$ "))
q <- seq(0.001, 9, by = 0.01)
lines(q, lnphden(q, rep(1,length(q) ), theta, alpha, S), col = "#3498DB", lwd = 2.5)


mu <- -0.35
sigma <- 0.8
thephfit <- ph(alpha = alpha, S = S)
p <- seq(0.001, 1-0.001, by = 0.001)
base::plot(quan(thephfit, p)$quantile, qlnorm(p, meanlog = mu, sdlog = sigma))
abline(0,1)

sequ <- seq(0.01, 6, by = 0.01)
base::plot(sequ, dlnorm(sequ, meanlog = mu, sdlog = sigma), type = "l")
lines(sequ, dens(thephfit, sequ)$dens, col= "red")

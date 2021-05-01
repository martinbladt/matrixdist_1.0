# Parameters of the matrix-Weibull
alpha0 <- c(0.5, 0.3, 0.2)
S0 <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
beta0 <- 2

ph <- ph(alpha = alpha0, S = S0)
iph <- iph(ph, gfun = "Weibull", gfun_pars = beta0)


sum(log(dens(iph,value)$dens ) * prob)
# Discretisation of the matrix-Weibull
truncationPoint <- 10
maxProbability <-  0.001
maxDeltat <- 0.05

deltat <- 0
t <- 0.0001 


prob = numeric(0)
value = numeric(0)

j <- 1

while (t < truncationPoint) {
  if (mWeibullden(t, alpha0, S0, beta0) < maxProbability / maxDeltat) {
    deltat = maxDeltat
  }
  else {
    deltat = maxProbability / mWeibullden(t, alpha0, S0, beta0)
  }
  proba_aux = deltat / 6 * (mWeibullden(t, alpha0, S0, beta0) + 4 * mWeibullden(t + deltat / 2, alpha0, S0, beta0) + mWeibullden(t + deltat, alpha0, S0, beta0))
  while (proba_aux > maxProbability) {
    deltat = deltat * 0.9
    proba_aux = deltat / 6 * (mWeibullden(t, alpha0, S0, beta0) + 4 * mWeibullden(t + deltat / 2, alpha0, S0, beta0) + mWeibullden(t + deltat, alpha0, S0, beta0))
  }
  if (proba_aux > 0) {
    value[j] = (t * mWeibullden(t, alpha0, S0, beta0) + 4 * (t + deltat / 2) * mWeibullden(t + deltat / 2, alpha0, S0, beta0) + (t + deltat) * mWeibullden(t + deltat, alpha0, S0, beta0)) / (mWeibullden(t, alpha0, S0, beta0) + 4 * mWeibullden(t + deltat / 2, alpha0, S0, beta0) + mWeibullden(t + deltat, alpha0, S0, beta0))
    prob[j] = proba_aux
    j <- j + 1
  }
  t = t + deltat
}

length(prob)
sum(prob)


# Function related to the fit
mix_dens <- function(z, pars){
  suppressWarnings(do.call(dpstable, append(pars,list(x = z))))
}
mLL <- function(parameters, valfn, L) {
  cum <- 0
  for (i in 1:(length(L))) {
    cum = cum - log(mix_dens(valfn[i], parameters)) * L[i]
  }
  return(cum)
}


#Initial values of the distribution be fitted
set.seed(1)
ph1 <- ph(structure = "GCoxian", dimension = 3)
pi <- ph1@pars$alpha
S <- ph1@pars$S
I <- diag(4)
alpha_o <- 1
kappa_o <- 1 
parameters <- c(alphaf,kappaf)
sum(log(dens(ph1, value)$dens) * prob)

density_isph_pstable <- function(x, alpha, S, par1, par2) {
  pars0 <- plogis(par1)
  den = numeric(0)
  for (i in 1:length(x)) {
    eg <- eigen(- x[i]^par2 * S)
    p <- eg$vectors
    d <- diag(exp(- (eg$values)^pars0) * (pars0 * (eg$values)^{pars0 - 1}) )
    mm <- p %*% d %*% solve(p)
    den[i] = sum(alpha %*%  mm %*% (S * (-1))) * par2 * x[i]^(par2-1)
  }
  return(den)
}

mLL2 <- function(x, weight, alpha, S, par1, par2) {
  -sum(log(density_isph_pstable(x, alpha, S, par1, par2)) * weight )
}

sum(density_isph_pstable(seq(0.01,10,by=0.01),pi, S, alpha_o, kappa_o) * 0.01)

sum(log(density_isph_pstable(value,pi, S, alpha_o, kappa_o)) * prob)
mLL2(value, prob, pi, S, alpha_o, kappa_o)

truncationPoint <- 10
maxProbability <-  0.01
maxDeltat <- 0.5
alpha_fit <- alpha_o
kappa_fit <- kappa_o

rcen = numeric(0)
rcenweight = numeric(0)
tl = numeric(0)

stepsEM <- 20
for (k in 1:stepsEM) {
  
  valuetrans <- value^(kappa_fit)
  
  theDensity <- density_cph_pstable(valuetrans, pi, S, alpha_fit)
  
  grid <- gridDensity_Rcpp(mix_dens, alpha_fit, truncationPoint, maxProbability, maxDeltat) 
  
  L <- EMstep_mPH_Pade(pi, S, valuetrans, prob,  theDensity, rcen, rcenweight, tl, grid$val, grid$prob)
  
  pi <- pi/sum(pi) #To compensate for the error due  to discretization - double check
  
  opt1 <- suppressWarnings(optim(par = alpha_fit, fn = mLL, valfn = grid$val, L = L, list(maxit = 50)))
  
  alpha_fit <- opt1$par
  
  opt2 <- suppressWarnings(optim(par = kappa_fit, fn = mLL2, x =value , weight = prob, alpha = pi, S = S, par1 = alpha_fit, list(maxit = 50)))
  kappa_fit <- opt2$par
  
  cat("\r", "iteration:", k,
      ", logLik:", - opt2$value,", factor:", kappa_fit * plogis(alpha_fit),
      sep = " ")
}

kappa_fit * plogis(alpha_fit)


sq <- seq(0.01, 8, by = 0.01)
base::plot(sq, dens(iph, sq)$dens , type = "l", xlab = "x", main = "Original density vs fitted density", ylab = "Density", lwd = 6,  cex.main=1.8, cex.lab=1.5, cex.axis=1.5, bty="l")
lines(sq, density_isph_pstable(sq, pi, S, alpha_fit, kappa_fit), lwd = 2, col = "#E74C3C")
legend("topright", 
       legend = c("Matrix-Weibull", 
                  "CPH with positive stable mixing"), 
       col = c("black", "#E74C3C"), 
       lty = c(2, 1), 
       bty = "n", 
       lwd = c(6,2), 
       cex = 1.1, 
       text.col = "black", 
       horiz = FALSE , 
       inset = c(0.05, 0.05))


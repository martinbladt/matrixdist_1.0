library(matrixdist)
library(survival)

veteran$status <- veteran$status^0

veteran <- rbind(veteran,veteran)

vfit <- coxph(Surv(time, status) ~ trt + karno + prior, veteran)
zp <- cox.zph(vfit, transform= function(time) log(time +20))

base::plot(zp[2]) # a plot for the 2nd variable in the fit
abline(0,0, col=2)
abline(h= vfit$coef[2], col=3, lwd=2, lty=2)

# time dependent coefficients can be used to address this..
# instead let's try our package

dat <- veteran
dat$trt <- as.numeric(dat$trt == 2)
dat$prior <- as.numeric(dat$prior == 10)
dat$karno <- dat$karno/10
dat$time <- dat$time/100

##

y <- dat$time[dat$status == 1]
rcen <- dat$time[dat$status == 0]
head(dat)
X <- dat[order(dat$status, decreasing = TRUE), c(1,5,8)] #8
head(dat)
set.seed(1)

dev.off()
A <- iph(ph(structure = "general", dimension = 2), gfun = "pareto", gfun_par = 1/2)

B <- fit(A, y = y, rcen = rcen, stepsEM = 1000, methods = c("RK", "RK"))
B <- fit(A, y = y, rcen = rcen, stepsEM = 1000, methods = c("RK", "UNI"))
B <- fit(A, y = y, rcen = rcen, stepsEM = 1000, methods = c("UNI", "UNI"), uni_epsilon = 0.0000000000000001)
B <- fit(A, y = y, rcen = rcen, stepsEM = 1000, methods = c("UNI", "PADE"))

base::plot(survival::survfit(survival::Surv(time, status) ~ 1, data = dat))
sq <- seq(0,10, by = .05)
pp <- matrixdist::cdf(B, sq, lower.tail = FALSE)
lines(sq, pp, col = "red")

################ now covariates

set.seed(1)
iA <- reg(x = iph(ph(structure = "coxian", dimension = 1),
                  gfun = "weibull", gfun_pars = 1),
          y = y, rcen = rcen, X = X, stepsEM = 300, methods = c("RK", "UNI"))

sqrt(diag(solve(-iA@fit$hessian)))
sqrt(diag(solve(length(y)*Fisher(x,y,X))))
-iA@fit$hessian[2:4,2:4]
Fisher(x,y,X)[1:3,1:3]
X

set.seed(1)
iB <- reg(x = iph(ph(structure = "general", dimension = 5),
                  gfun = "weibull", gfun_pars = 1),
          y = y, rcen = rcen, X = X, stepsEM = 300,
          methods = c("RK","UNI"))

set.seed(1)
iB <- reg(x = iph(ph(structure = "coxian", dimension = 3),
                  gfun = "weibull", gfun_pars = 1),
          y = y, rcen = rcen, X = X, stepsEM = 300,
          methods = c("RK","PADE"))
set.seed(1)
iB <- reg(x = iph(ph(structure = "general", dimension = 5),
                  gfun = "weibull", gfun_pars = 1),
          y = y, rcen = rcen, X = X, stepsEM = 1000,
          methods = c("UNI","UNI"), uni_epsilon = 0.0000000000000001)

stats4::AIC(iB)

summary(survival::survreg(survival::Surv(time,status) ~  trt + karno+prior ,dist="exponential", data = dat))
logLik(survival::survreg(survival::Surv(time,status) ~  trt + karno ,dist="exponential", data = dat))
logLik(survival::survreg(survival::Surv(time,status) ~  trt + karno ,dist="weibull", data = dat))


#### A PLOT
par(mfrow=c(3,3))
for(i in 2:9){
  x0 <- c(0, i, 1)
  dat_sub <- subset(dat, trt == x0[1] & karno == x0[2])
  iAplot <- evaluate(iA, subject = x0)
  iBplot <- evaluate(iB, subject = x0)
  obj <- survival::survfit(survival::Surv(time,status) ~  1, data = dat_sub)
  
  base::plot(log(obj$time),obj$surv,
             xlab = "time", 
             ylab = "Survival Probability",
             ylim = c(0,1))
  sq <- seq(0,max(obj$time), by = .001)
  lines(log(sq), (matrixdist::cdf(iAplot, sq, lower.tail = FALSE)), col = "red")
  lines(log(sq), (matrixdist::cdf(iBplot, sq, lower.tail = FALSE)), col = "blue")
}

## cox-snell residuals comparision
times <- c(y, rcen)

res1 <- numeric(nrow(dat))
res2 <- numeric(nrow(dat))

for(i in 1:nrow(dat)){
  res1[i] <- -log(matrixdist::cdf(evaluate(iA, subject = unlist(X[i,])), times[i], lower.tail = FALSE))
  res2[i] <- -log(matrixdist::cdf(evaluate(iB, subject = unlist(X[i,])), times[i], lower.tail = FALSE))
}
dat$coxsnell1 <- res1
dat$coxsnell2 <- res2

## GOF
labels = c("a", "b")
par(mfrow=c(2,1))
for(i in 9:10){
  m1 <- survfit(survival::Surv(dat[,i],dat$status) ~  1)
  base::plot(log(m1$time), m1$surv, type = "s", ylab = "Survival Probability", xlab = "log-time",
             main = paste("Model", labels[i-8]))
  lines(log(m1$time), m1$upper, type = "s", lty = 3)
  lines(log(m1$time), m1$lower, type = "s", lty = 3)
  sq <- exp(seq(-10,3,0.1))
  lines(log(sq), exp(-sq), col = "red", lty = 2)
}

##
## Fit model on Cox-Snell residuals (Approximately Expo(1) distributed under correct model)

## Nelson-Aalen estimator for baseline hazard (all covariates zero)
par(mfrow=c(2,1))
for(i in 9:10){
  df_base_haz <- basehaz(coxph(formula = Surv(dat[,i],dat$status) ~ 1),centered = FALSE)
  base::plot((df_base_haz$time), (df_base_haz$hazard),
             xlab = "time", ylab = "Cumulative Hazard",
             main = paste("Model", labels[i-8]))
  abline(0,1)
}



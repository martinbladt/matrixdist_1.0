library(matrixdist)
library(survival)

dim(veteran)

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

set.seed(1)

A <- iph(ph(structure = "Coxian", dimension = 2), gfun = "Pareto", gfun_par = 1/2)

B <- fit(A, y = y, rcen = rcen, stepsEM = 2000)

base::plot(survival::survfit(survival::Surv(time, status) ~ 1, data = dat))
sq <- seq(0,10, by = .05)
pp <- matrixdist::cdf(B, sq, lower.tail = FALSE)
lines(pp$q, pp$cdf, col = "red")

################ now covariates

set.seed(1)
iA <- reg(x = iph(ph(structure = "Coxian", dimension = 1),
                  gfun = "Weibull", gfun_pars = 1),
          y = y, rcen = rcen, X = X, stepsEM = 300)

set.seed(2)
iB <- aft(x = iph(ph(structure = "Coxian", dimension = 1), 
                  gfun = "LogLogistic", gfun_pars = c(1,1)),
          y = y, rcen = rcen, X = X, stepsEM = 500)

set.seed(1)
iC <- reg(x = iph(ph(structure = "Coxian", dimension = 2), 
                  gfun = "Weibull", gfun_pars = 1),
          y = y, rcen = rcen, X = X, stepsEM = 1000)

set.seed(1)
iD <- aft(x = iph(ph(structure = "Coxian", dimension = 2), 
                  gfun = "LogNormal", gfun_pars = c(1)),
          y = y, rcen = rcen, X = X, stepsEM = 1000)

iA@fit$loglik
iB@fit$loglik
iC@fit$loglik
iD@fit$loglik

coef(iA)
coef(iB)
coef(iC)
coef(iD)

logLik(survival::survreg(survival::Surv(time,status) ~  trt + karno ,dist="exponential", data = dat))
logLik(survival::survreg(survival::Surv(time,status) ~  trt + karno ,dist="weibull", data = dat))


#### A PLOT
par(mfrow=c(3,3))
for(i in 2:9){
  x0 <- c(0, i)
  dat_sub <- subset(dat, trt == x0[1] & karno == x0[2])
  
  iAplot <- eval(iA, subject = x0)
  iBplot <- eval(iB, subject = x0)
  iCplot <- eval(iC, subject = x0)
  iDplot <- eval(iD, subject = x0)
  
  obj <- survival::survfit(survival::Surv(time,status) ~  1, data = dat_sub)
  
  base::plot(log(obj$time),obj$surv,
             xlab = "time", 
             ylab = "Survival Probability",
             ylim = c(0,1))
  sq <- seq(0,max(obj$time), by = .001)
  lines(log(sq), (matrixdist::cdf(iAplot, sq, lower.tail = FALSE)$cdf), col = "red")
  lines(log(sq), (matrixdist::cdf(iBplot, sq, lower.tail = FALSE)$cdf), col = "blue")
  lines(log(sq), (matrixdist::cdf(iCplot, sq, lower.tail = FALSE)$cdf), col = "orange")
  lines(log(sq), (matrixdist::cdf(iDplot, sq, lower.tail = FALSE)$cdf), col = "green")
}

# another plot

par(mfrow=c(1,2))
for(i in 0:1){
  x0 <- c(i, 9)
  dat_sub <- subset(dat, trt == x0[1] & karno == x0[2])
  iAplot <- eval(iA, subject = x0)
  iBplot <- eval(iB, subject = x0)
  iCplot <- eval(iC, subject = x0)
  iDplot <- eval(iD, subject = x0)
  
  obj <- survival::survfit(survival::Surv(time,status) ~  1, data = dat_sub)
  
  base::plot(log(obj$time),obj$surv,
             xlab = "time", 
             ylab = "Survival Probability",
             ylim = c(0,1))
  sq <- seq(0,max(obj$time), by = .001)
  lines(log(sq), (matrixdist::cdf(iAplot, sq, lower.tail = FALSE)$cdf), col = "red")
  lines(log(sq), (matrixdist::cdf(iBplot, sq, lower.tail = FALSE)$cdf), col = "blue")
  lines(log(sq), (matrixdist::cdf(iCplot, sq, lower.tail = FALSE)$cdf), col = "orange")
  lines(log(sq), (matrixdist::cdf(iDplot, sq, lower.tail = FALSE)$cdf), col = "green")
}

## cox-snell residuals comparision
times <- c(y, rcen)

res1 <- numeric(nrow(dat))
res2 <- numeric(nrow(dat))
res3 <- numeric(nrow(dat))
res4 <- numeric(nrow(dat))


for(i in 1:nrow(dat)){
  res1[i] <- -log(matrixdist::cdf(eval(iA, subject = unlist(X[i,])), times[i], lower.tail = FALSE)$cdf)
  res2[i] <- -log(matrixdist::cdf(eval(iB, subject = unlist(X[i,])), times[i], lower.tail = FALSE)$cdf)
  res3[i] <- -log(matrixdist::cdf(eval(iC, subject = unlist(X[i,])), times[i], lower.tail = FALSE)$cdf)
  res4[i] <- -log(matrixdist::cdf(eval(iD, subject = unlist(X[i,])), times[i], lower.tail = FALSE)$cdf)
}
dat$coxsnell1 <- res1
dat$coxsnell2 <- res2
dat$coxsnell3 <- res3
dat$coxsnell4 <- res4

## GOF
labels = c("a", "b", "c", "d")
png("/Users/martinbladt/Dropbox/Regression for PH/gof_vet_surv.png",
    width = 12, height = 12, units = 'in', res = 200)
par(mfrow=c(2,2))
for(i in 9:12){
  m1 <- survfit(survival::Surv(dat[,i],dat$status) ~  1)
  base::plot(log(m1$time), m1$surv, type = "s", ylab = "Survival Probability", xlab = "log-time",
             main = paste("Model", labels[i-8]))
  lines(log(m1$time), m1$upper, type = "s", lty = 3)
  lines(log(m1$time), m1$lower, type = "s", lty = 3)
  sq <- exp(seq(-10,3,0.1))
  lines(log(sq), exp(-sq), col = "red", lty = 2)
}
dev.off()

##
## Fit model on Cox-Snell residuals (Approximately Expo(1) distributed under correct model)

## Nelson-Aalen estimator for baseline hazard (all covariates zero)
png("/Users/martinbladt/Dropbox/Regression for PH/gof_vet_haz.png",
    width = 12, height = 12, units = 'in', res = 200)
par(mfrow=c(2,2))
for(i in 9:12){
  df_base_haz <- basehaz(coxph(formula = Surv(dat[,i],dat$status) ~ 1),centered = FALSE)
  base::plot((df_base_haz$time), (df_base_haz$hazard),
             xlab = "time", ylab = "Cumulative Hazard",
             main = paste("Model", labels[i-8]))
  abline(0,1)
}
dev.off()



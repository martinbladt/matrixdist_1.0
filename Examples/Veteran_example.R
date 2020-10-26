library(matrixdist)
library(survival)

dim(veteran)

vfit <- coxph(Surv(time, status) ~ trt + prior + karno, veteran)
zp <- cox.zph(vfit, transform= function(time) log(time +20))

base::plot(zp[3]) # a plot for the 3rd variable in the fit
abline(0,0, col=2)
abline(h= vfit$coef[3], col=3, lwd=2, lty=2)

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

X <- dat[order(dat$status, decreasing = TRUE), c(1,5,8)]

set.seed(1)

A <- ph(structure = "Coxian", dimension = 2)

B <- fit(A, y = y, rcen = rcen, stepsEM = 2000)

base::plot(survival::survfit(survival::Surv(time, status) ~ 1, data = dat))
sq <- seq(0,10, by = .05)
pp <- matrixdist::cdf(B, sq, lower.tail = FALSE)
lines(pp$q, pp$cdf, col = "red")

################ now covariates

set.seed(1)
A <- ph(structure = "Coxian", dimension = 1)
iA <- iph(A, gfun = "Weibull", gfun_pars = 1) 
iB <- reg(x = iA, y = y, rcen = rcen, X = X, stepsEM = 500)
#iC <- aft(x = iA, y = y, rcen = rcen, X = X, stepsEM = 500)
iB@fit

set.seed(1)
A <- ph(structure = "Coxian", dimension = 2)
iA <- iph(A, gfun = "Weibull", gfun_pars = 1) 
iC <- reg(x = iA, y = y, rcen = rcen, X = X, stepsEM = 500)
iC@fit

#iD <- reg2(x = iA, y = y, rcen = rcen, X = X, stepsEM = 500)

logLik(survival::survreg(survival::Surv(time,status) ~  trt + prior + karno ,dist="exponential", data = dat))
logLik(survival::survreg(survival::Surv(time,status) ~  trt + prior + karno ,dist="weibull", data = dat))


coef(iB)
coef(iC)
#coef(iD)

#### A PLOT
x0 <- c(1, 7, 0)
dat_sub <- subset(dat, trt == x0[1] & karno == x0[2] & prior == x0[3])

iBplot <- eval(iB, subject = x0)
iCplot <- eval(iC, subject = x0)
#iDplot <- eval(iD, subject = x0)


obj <- survival::survfit(survival::Surv(time,status) ~  1, data = dat_sub)

base::plot(log(obj$time),obj$surv,
           xlab = "time", 
           ylab = "Cumulative Hazard",
           ylim = c(0,1))
sq <- seq(0,max(obj$time), by = .001)
lines(log(sq), (matrixdist::cdf(iBplot, sq, lower.tail = FALSE)$cdf), col = "red")
lines(log(sq), (matrixdist::cdf(iCplot, sq, lower.tail = FALSE)$cdf), col = "orange")
#lines(log(sq), (cdf(iDplot, sq, lower.tail = FALSE)$cdf), col = "blue")

## cox-snell residuals comparision
times <- c(y, rcen)

res1 <- numeric(nrow(dat))
res2 <- numeric(nrow(dat))

for(i in 1:nrow(dat)){
  res1[i] <- -log(matrixdist::cdf(eval(iB, subject = unlist(X[i,])), times[i], lower.tail = FALSE)$cdf)
  res2[i] <- -log(matrixdist::cdf(eval(iC, subject = unlist(X[i,])), times[i], lower.tail = FALSE)$cdf)
}
dat$coxsnell1 <- res1
dat$coxsnell2 <- res2

##
## Fit model on Cox-Snell residuals (Approximately Expo(1) distributed under correct model)
fit_coxsnell1 <- coxph(formula = Surv(coxsnell1, status) ~ 1, data    = dat)
fit_coxsnell2 <- coxph(formula = Surv(coxsnell2, status) ~ 1, data    = dat)

## Nelson-Aalen estimator for baseline hazard (all covariates zero)
df_base_haz1 <- basehaz(fit_coxsnell1, centered = FALSE)
df_base_haz2 <- basehaz(fit_coxsnell2, centered = FALSE)

base::plot(df_base_haz1)
abline(0,1)

base::plot(df_base_haz2)
abline(0,1)




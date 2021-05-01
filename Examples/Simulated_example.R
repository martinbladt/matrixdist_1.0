# Simulation example
library(matrixdist)
n0 <- 500

dat <- data.frame(sex = c(rep(1, n0), rep(0, n0)))
A1 <- ph(alpha = c(1/4, 1/2, 1/4), 
         S = matrix(c(-10, 0, 0 ,0 , -1, 0, 0, 0, -1/10), 3, 3))
iA1 <- iph(A1,gfun = "Weibull", gfun_pars = 1.5)
A2 <- A1; A2@pars$S <- 2 * A2@pars$S
iA2 <- iph(A2,gfun = "Weibull", gfun_pars = 1)

set.seed(1000) #500 is good
group1 <- sim(iA1, n0)
group2 <- sim(iA2, n0)

censoring <- rexp(2*n0, rate = 1/10) #1/5
dat$time <- pmin(c(group1, group2), censoring)
dat$status <- as.numeric(dat$time != censoring)

head(dat)

y <- dat$time[dat$status == 1]
rcen <- dat$time[dat$status == 0]

X <- dat[order(dat$status, decreasing = TRUE), c(1)]

### fitting regressions

set.seed(1)
A <- ph(structure = "Coxian", dimension = 1)
iA <- iph(A, gfun = "Weibull", gfun_pars = 1) 
iB <- reg(x = iA, y = y, rcen = rcen, X = X, stepsEM = 100)

set.seed(1)
A <- ph(structure = "Coxian", dimension = 1)
iA <- iph(A, gfun = "Weibull", gfun_pars = 1) 
iC <- reg2(x = iA, y = y, rcen = rcen, X = X, stepsEM = 100)

iD <- reg2(x = iA2, y = y, rcen = rcen, X = X, stepsEM = 150)

#### A PLOT
x0 <- c(1)
dat_sub <- subset(dat, sex == x0)

iBplot <- eval(iB, subject = x0)
iCplot <- eval(iC, subject = x0)
iDplot <- eval(iD, subject = x0)

iDplot@gfun$pars

obj <- survival::survfit(survival::Surv(time,status) ~  1, data = dat_sub)

#png("/Users/martinbladt/Dropbox/Regression for PH/sim_g2.png",
#    width = 6, height = 6, units = 'in', res = 200)
par(mfrow=c(1,1))
base::plot(log(obj$time), obj$surv,
           xlab = "log-time", 
           ylab = "Survival Probability",
           main = "Group 2",
           ylim = c(0,1), pch = 20)
sq <- exp(seq(-50, 2* log(max(obj$time)), by = .01))
lines(log(sq), (matrixdist::cdf(iBplot, sq, lower.tail = FALSE)$cdf), col = "red", lwd = 2, lty = 2)
lines(log(sq), (matrixdist::cdf(iCplot, sq, lower.tail = FALSE)$cdf), col = "orange", lwd = 2, lty = 3)
lines(log(sq), (matrixdist::cdf(iDplot, sq, lower.tail = FALSE)$cdf), col = "purple", lwd =2)
dev.off()
## estimated vs theoretical hazard ratio


#png("/Users/martinbladt/Dropbox/Regression for PH/hazard_ratio.png",
#width = 6, height = 6, units = 'in', res = 200)
par(mfrow=c(1,1))
sq <- seq(0.01, 10, by = 0.01)
ratio <- haz(iA1, sq)$haz/haz(iA2, sq)$haz
base::plot(sq, ratio, type = "l", lty = 1, lwd = 2, main = "Hazard ratio", ylab = "Ratio", xlab = "time")
lines(sq, haz(eval(iB, subject = 1), sq)$haz/haz(eval(iB, subject = 0), sq)$haz, lty = 2, col = "red", lwd = 1.5)
lines(sq, haz(eval(iC, subject = 1), sq)$haz/haz(eval(iC, subject = 0), sq)$haz, lty = 3, col = "orange", lwd = 1.5)
lines(sq, haz(eval(iD, subject = 1), sq)$haz/haz(eval(iD, subject = 0), sq)$haz, lty = 1, col = "purple", lwd = 1.5)
dev.off()

## cox-snell residuals comparision
times <- c(y, rcen)

res1 <- numeric(nrow(dat))
res2 <- numeric(nrow(dat))
res3 <- numeric(nrow(dat))

for(i in 1:nrow(dat)){
  res1[i] <- -log(matrixdist::cdf(matrixdist::eval(iB, subject = unlist(X[i])), times[i], lower.tail = FALSE)$cdf)
  res2[i] <- -log(matrixdist::cdf(matrixdist::eval(iC, subject = unlist(X[i])), times[i], lower.tail = FALSE)$cdf)
  res3[i] <- -log(matrixdist::cdf(matrixdist::eval(iD, subject = unlist(X[i])), times[i], lower.tail = FALSE)$cdf)
}
dat$coxsnell1 <- res1
dat$coxsnell2 <- res2
dat$coxsnell3 <- res3

##
## Fit model on Cox-Snell residuals (Approximately Expo(1) distributed under correct model)
fit_coxsnell1 <- coxph(formula = Surv(coxsnell1, status) ~ 1, data = dat)
fit_coxsnell2 <- coxph(formula = Surv(coxsnell2, status) ~ 1, data = dat)
fit_coxsnell3 <- coxph(formula = Surv(coxsnell3, status) ~ 1, data = dat)

## Nelson-Aalen estimator for baseline hazard (all covariates zero)
df_base_haz1 <- basehaz(fit_coxsnell1, centered = FALSE)
df_base_haz2 <- basehaz(fit_coxsnell2, centered = FALSE)
df_base_haz3 <- basehaz(fit_coxsnell3, centered = FALSE)



#png("/Users/martinbladt/Dropbox/Regression for PH/gof_QQ.png",
#    width = 18, height = 6, units = 'in', res = 200)
par(mfrow=c(1,3))
base::plot(log(df_base_haz1$time), log(df_base_haz1$hazard),
           xlab = "log-time", ylab = "log-cumhaz")
abline(0,1)
base::plot(log(df_base_haz2$time), log(df_base_haz2$hazard),
           xlab = "log-time", ylab = "log-cumhaz")
abline(0,1)
base::plot(log(df_base_haz3$time), log(df_base_haz3$hazard),
           xlab = "log-time", ylab = "log-cumhaz")
abline(0,1)
dev.off()


#png("/Users/martinbladt/Dropbox/Regression for PH/gof_surv.png",
#    width = 18, height = 6, units = 'in', res = 200)
par(mfrow=c(1,3))
m1 <- survfit(survival::Surv(coxsnell1,status) ~  1, data = dat)
base::plot(log(m1$time), m1$surv, type = "s", ylab = "Survival Probability", xlab = "log-time")
lines(log(m1$time), m1$upper, type = "s", lty = 3)
lines(log(m1$time), m1$lower, type = "s", lty = 3)
sq <- exp(seq(-10,3,0.1))
lines(log(sq), exp(-sq), col = "red", lty = 2)
m2 <- survfit(survival::Surv(coxsnell2,status) ~  1, data = dat)
base::plot(log(m2$time), m2$surv, type = "s", ylab = "Survival Probability", xlab = "log-time")
lines(log(m2$time), m2$upper, type = "s", lty = 3)
lines(log(m2$time), m2$lower, type = "s", lty = 3)
sq <- exp(seq(-10,3,0.1))
lines(log(sq), exp(-sq), col = "red", lty = 2)
m3 <- survfit(survival::Surv(coxsnell3,status) ~  1, data = dat)
base::plot(log(m3$time), m3$surv, type = "s", ylab = "Survival Probability", xlab = "log-time")
lines(log(m3$time), m3$upper, type = "s", lty = 3)
lines(log(m3$time), m3$lower, type = "s", lty = 3)
sq <- exp(seq(-10,3,0.1))
lines(log(sq), exp(-sq), col = "red", lty = 2)
dev.off()
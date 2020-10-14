data("lung", package = "survival")

standardize <- function(x) x/max(x)

lung$status <- as.numeric(lung$status > 1)
lung$sex <- as.numeric(lung$sex > 1)
lung$age <- standardize(lung$age)
lung$time <- standardize(lung$time)



y <- lung$time[lung$status == 1]
rcen <- lung$time[lung$status == 0]

X <- lung[order(lung$status, decreasing = TRUE), c(5)] #(4,5) is good


set.seed(1)

A <- ph(structure = "Coxian", dimension = 4)

B <- fit(A, y = y, rcen = rcen, stepsEM = 2000)

plot(survival::survfit(survival::Surv(time, status) ~ 1, data = lung), 
     xlab = "Days", 
     ylab = "Overall survival probability")
sq <- seq(0,1, by = .05)
pp <- cdf(B, sq, lower.tail = FALSE)
lines(pp$q, pp$cdf, col = "red")

################ now covariates

set.seed(1)
A <- ph(structure = "Coxian", dimension = 3)

#B <- reg(x = A, y = y, rcen = rcen, X = X)
#C <- aft(x = A, y = y, rcen = rcen, X = X) #gives the same, as it should

#logLik(survival::survreg(survival::Surv(time,status) ~  age + sex,dist="exponential", data = lung))

iA <- iph(A, gfun = "Weibull", gfun_pars = 1) 
iB <- reg(x = iA, y = y, rcen = rcen, X = X, stepsEM = 500)

iD <- reg2(x = iA, y = y, rcen = rcen, X = X, stepsEM = 500)

#logLik(survival::survreg(survival::Surv(time,status) ~  age + sex,dist="weibull", data = lung))

#### A PLOT
x0 <- c(0)
lung_sub <- subset(lung, sex == x0)

iBplot <- iB; iBplot@pars$S <- iBplot@pars$S * exp(sum(c(-0.5621967 )*x0))
iDplot <- iD; iDplot@pars$S <- iDplot@pars$S * exp(sum(c(-0.2304824)*x0)) 
iDplot@gfun$pars <- exp(0.1336482 + sum(c(0.2357306)*x0))

obj <- survival::survfit(survival::Surv(time,status) ~  1, data = lung_sub)

plot(log(obj$time),obj$surv,
     xlab = "Days", 
     ylab = "Cumulative Hazard",
     ylim = c(0,1))
sq <- seq(0,max(obj$time), by = .001)
lines(log(sq), (cdf(iBplot, sq, lower.tail = FALSE)$cdf), col = "red")
lines(log(sq), (cdf(iDplot, sq, lower.tail = FALSE)$cdf), col = "blue")
####

res.cox <- survival::coxph(survival::Surv(time, status) ~ ., data = lung)
survival::cox.zph(res.cox)

head(lung)

###

iA <- iph(A, gfun = "Pareto", gfun_pars = 1) 

iB <- reg(x = iA, y = y, rcen = rcen, X = X)
iC <- aft(x = iA, y = y, rcen = rcen, X = X)

iA <- iph(A, gfun = "Gompertz", gfun_pars = 1) 

iB <- reg(x = iA, y = y, rcen = rcen, X = X)
iC <- aft(x = iA, y = y, rcen = rcen, X = X)


iA <- iph(A, gfun = "LogLogistic", gfun_pars = c(1, 1)) 

iB <- reg(x = iA, y = y, rcen = rcen, X = X)
iB <- aft(x = iA, y = y, rcen = rcen, X = X)

# make regression objects and methods, optimize code everywhere, polish presentation and defaults to fit



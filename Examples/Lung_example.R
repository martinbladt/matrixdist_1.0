data("lung", package = "survival")

standardize <- function(x) x/max(x)

lung$status <- as.numeric(lung$status > 1)
lung$sex <- as.numeric(lung$sex > 1)
lung$age <- standardize(lung$age)
lung$time <- standardize(lung$time)



y <- lung$time[lung$status == 1]
rcen <- lung$time[lung$status == 0]

X <- lung[order(lung$status, decreasing = TRUE), c(4, 5)]


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
A <- ph(structure = "Coxian", dimension = 2)

B <- reg(x = A, y = y, rcen = rcen, X = X)
C <- aft(x = A, y = y, rcen = rcen, X = X) #gives the same, as it should

iA <- iph(A, gfun = "Weibull", gfun_pars = 1) 

iB <- reg(x = iA, y = y, rcen = rcen, X = X)
iC <- aft(x = iA, y = y, rcen = rcen, X = X) #gives the same, as it should

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



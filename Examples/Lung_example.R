data("lung", package = "survival")
head(lung)

plot(lung$time)

set.seed(1)

A <- ph(structure = "Coxian", dimension = 2)

#B <- fit(A, lung$time, stepsEM = 2000)

#m_plot(B, lung$time)

y <- lung$time[lung$status == 2]
rcen <- lung$time[lung$status == 1]

B <- fit(A, y = y, rcen = rcen, stepsEM = 2000)

plot(survival::survfit(survival::Surv(time, status) ~ 1, data = lung), 
     xlab = "Days", 
     ylab = "Overall survival probability")
sq <- seq(0,1100, by = 1)
pp <- cdf(B, sq, lower.tail = FALSE)
lines(pp$q, pp$cdf, col = "red")



################ now covariates
X <- rev(lung[order(lung$status), c(4, 5)])/100 #clean this up
X <- t(t(X)-colMeans(X)) #clean this up

set.seed(10)
A <- ph(structure = "Coxian", dimension = 2)

B <- reg(x = A, y = y, rcen = rcen, X = X, stepsEM = 1000)
m_plot(B, y)

iA <- iph(A, gfun = "Pareto", gfun_pars = 1) 

iB <- reg(x = iA, y = y, rcen = rcen, X = X, stepsEM = 1000) #no jala!
m_plot(B, y)



# make regression objects and methods














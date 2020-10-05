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
X <- rev(lung[order(lung$status), c(4, 5)])

B <- reg(x = A, y = y/100, rcen = rcen/100, X = X/100, stepsEM = 1000)

m_plot(B, y/100)

B #degenerated to zero! likelihood increased and then decreased :(


head(lung)














set.seed(123456)
ph1 <- ph(structure = "Coxian", dimension = 3)
cph1 <- cph(ph1, dens = "gamma", dens_pars = c(shape = 2, rate = 1)) #first parameter is alpha
#cph1 <- cph(ph1, dens = "pstable", dens_pars = c(frac = 3))
data0 <- sim(cph1, 100)
data <- pmin(data0, 4)

y <- data[data < 4]
rcen <- data[data == 4]


cph1
coef(cph1)

dens(cph1, 1:10)$dens
cdf(cph1, 1:10)
cdf(cph1, 1:10, lower.tail = F)
quan(cph1, c(0.5, 0.6, 0.7))
haz(cph1, 1:10)

f <- fit(x = cph1, y = data, rcen = rcen, stepsEM = 20)

#cph2 <- cph(ph1, dens = "gamma", dens_pars = c(shape = 2)) #first parameter is alpha
#f2 <- fit(x = cph2, y = data, stepsEM = 200) 

evir::hill(data, start = 5, option = "xi")

hist(log(data), breaks = 50, freq = FALSE)
sq <- seq(-20, 10, by = 0.1)
lines(sq, matrixdist::dens(cph1, exp(sq))$dens*exp(sq), lwd = 2)
lines(sq, matrixdist::dens(f, exp(sq))$dens*exp(sq), col = "red", lwd = 2)
#lines(sq, matrixdist::dens(f2, exp(sq))$dens*exp(sq), col = "green", lwd = 2)


set.seed(123456)
ph1 <- ph(structure = "Coxian", dimension = 1)

cph1 <- cph(ph1, dens = "gamma", dens_pars = c(shape = 2, rate = 1)) #first parameter is alpha

cph1 <- cph(ph1, dens = "norm", dens_pars = c(mean = 0, sd = 1)) 

cph1 <- cph(ph1, dens = "unif", dens_pars = c(max = 4)) 

cph1 <- cph(ph1, dens = "beta", dens_pars = c(shape1 = 2, shape2 = 2)) #first parameter is alpha

cph1 <- cph(ph1, dens = "lnorm", dens_pars = c(meanlog = 1, sdlog = 8)) #first parameter is alpha

cph1 <- cph(ph1, dens = "weibull", dens_pars = c(shape = 0.1, scale = 1)) #first parameter is alpha

data <- sim(cph1, 2000)

base::plot(data)

evir::hill(data, start = 5, option = "xi")
evir::meplot(data)

# weights and censoring both missing...

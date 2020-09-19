# An example

alpha <- c(0.5, 0.3, 0.2)
S <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)

A <- ph(alpha, S)
d(A, y = c(1, 2, 3))
p(A, q = c(1, 2, 3))


A

data <- r(A)

m_plot(A, data)
m_plot(A)

B <- fit(A, data)

B
m_plot(B, data)
m_plot(B)

# Inhomogenous

A <- ph(structure = "Coxian")

iA <- iph(A, gfun = "Weibull", gfun_pars = 2)

iA

data <- r(iA)

plot(data)

evmix::hillplot(data)

iB <- fit(iA, data, stepsEM = 400)

iB

iC <- iph(A, gfun = "Pareto", gfun_pars = 2)

iC

data <- r(iC)

plot(data)

evmix::hillplot(data)

iD <- fit(iC, data, stepsEM = 300)

iD

iE <- iph(A, gfun = "Gompertz", gfun_pars = 2)

iE

data <- r(iE)

plot(data)

evmix::hillplot(data)

iF <- fit(iE, data, stepsEM = 300)

iF

iG <- iph(A, gfun = "GEVD", gfun_pars = c(1, 1, 2)) #this crashes R

iG

data <- r(iG)

plot(data)

evmix::hillplot(data)

x <- iG;y <- data
iH <- fit(iG, data, stepsEM = 300)

iH

## clone par_g??

# MPH*

alpha <- c(0.5, 0.3, 0.2)
S <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
R <- matrix(cbind(c(1,0,0.8),c(0,1,0.2)), nrow = 3, ncol = 2)

A <- ph(alpha, S)
A

B <- mph(A, R)
B

r(B, n = 10)

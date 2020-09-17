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

iB <- fit(iA, data, stepsEM = 100)

iB

# MPH*

alpha <- c(0.5, 0.3, 0.2)
S <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
R <- matrix(cbind(c(1,0,0.8),c(0,1,0.2)), nrow = 3, ncol = 2)

A <- ph(alpha, S)
A

B <- mph(A, R)
B

r(B, n = 10)

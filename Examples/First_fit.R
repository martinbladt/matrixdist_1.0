# An example

alpha <- c(0.5, 0.3, 0.2)
S <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)

A <- ph(alpha, S)
d(A, y = c(1, 2, 3))

A

data <- r(A)

m_plot(A, data)
m_plot(A)

B <- fit(A, data)

B
m_plot(B, data)
m_plot(B)

# Inhomogenous

A <- ph(structure = "General")

iA <- iph(A, gfun = "Pareto", gfun_pars = 4)

iA

data <- r(iA)

plot(data)

evmix::hillplot(data)

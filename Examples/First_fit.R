# An example

alpha <- c(0.5, 0.3, 0.2)
S <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)

A <- ph(alpha, S)
d(A, y = c(1, 2, 3))

A

data <- r(A)

m_plot(A, data)
m_plot(A)

A_copy <- A

B <- fit(A, data)

B
m_plot(B, data)
m_plot(B) # very good


A #but this changed too..
A_copy # and even this, that didn't enter the function. the rcpp code seems to create pointers that spread!






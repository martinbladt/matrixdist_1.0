pi1 <- c(1,0)
S1 <- matrix(c(-1,0, 0.5, -0.5), nrow=2)

pi2 <- c(0.5, 0.5)
S2 <- matrix(c(-0.1,0, 0, -0.2), nrow=2)

sumPH(pi1, S1, pi2, S2)

sumPH(pi2, S2, pi1, S1)


A <- matrix(c(1,3, 2, 4), nrow=2)
B <- matrix(c(0,6, 5, 7), nrow=2)

Kroneckerproduct(A,B)
kronecker(A,B) #Ya estaba en R


A <- matrix(c(1,3, 2, 4), nrow=2)
B <- matrix(c(5,7, 10, 9), nrow=2)
Kroneckersum(A,B)

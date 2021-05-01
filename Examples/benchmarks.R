
library(microbenchmark)
set.seed(123)
M1 <- matrix(sample(1e4),ncol=100)
M2 <- matrix(sample(1e4),nrow=100)

identical(productArma(M1,M2), matrix_product(M1,M2) , M1 %*% M2)

microbenchmark(productArma(M1,M2),  matrix_product(M1,M2), M1 %*% M2, times=10000L)

identical(sumArma(M1,M2),  matrix_sum(M1,M2), M1 + M2)
microbenchmark(sumArma(M1,M2),  matrix_sum(M1,M2), M1 + M2, times=10000L)


A <- random_structure(100)$T

identical(invArma(A), matrix_inverse(A))

microbenchmark(invArma(A), matrix_inverse(A), times=10000L)

library(expm)
matrix_exponential(A)
mexponentialArma(A)
expm(A)

identical(round(mexponentialArma(A), digits = 10), round(matrix_exponential(A), digits = 10), round(expm(A), digits = 10))
microbenchmark(mexponentialArma(A), matrix_exponential(A), expm(A), times=1000L)




thePH <- ph(structure = "General", dimension = 4)

n <- 20000
set.seed(1)
x <- rphasetype(n, coef(thePH)$alpha, coef(thePH)$S/200) 
sum(log(phdensity(x, coef(thePH)$alpha, coef(thePH)$S/200)))
sum(log(phdensityArma(x, coef(thePH)$alpha, coef(thePH)$S/200)))


x <- seq(0.01,10, by = 0.001)
sum(log(phdensity(x, coef(thePH)$alpha, coef(thePH)$S)))
sum(log(phdensityArma(x, coef(thePH)$alpha, coef(thePH)$S)))

identical(round(phdensity(x, coef(thePH)$alpha, coef(thePH)$S), digits = 10), round(phdensityArma(x, coef(thePH)$alpha, coef(thePH)$S), digits = 10))
microbenchmark(phdensity(x, coef(thePH)$alpha, coef(thePH)$S), phdensityArma(x, coef(thePH)$alpha, coef(thePH)$S), times=1000L)


pars<-c(0.5, 1)
theDens <- function(z, pars){
  do.call(dgamma, append(pars,list(x = z)))
}

frac <- 0.8
gridDensity(dpstable, frac, truncationPoint = 10, maxProbability = 0.01, maxDeltat = 0.5)
gridDensity_Rcpp(dpstable, frac, truncationPoint = 10, maxProbability = 0.01, maxDeltat = 0.5)

identical(round(gridDensity(dpstable, frac, truncationPoint = 10, maxProbability = 0.01, maxDeltat = 0.5), digits = 10)
          ,round(gridDensity_Rcpp(dpstable, frac, truncationPoint = 10, maxProbability = 0.01, maxDeltat = 0.5), digits = 10))

microbenchmark(gridDensity(dpstable, frac, truncationPoint = 10, maxProbability = 0.01, maxDeltat = 0.5)
               ,gridDensity_Rcpp(dpstable, frac, truncationPoint = 10, maxProbability = 0.01, maxDeltat = 0.5), times=100L )

gridDensity(theDens, pars, truncationPoint = 10, maxProbability = 0.01, maxDeltat = 0.5)
gridDensity_Rcpp(theDens, pars, truncationPoint = 10, maxProbability = 0.01, maxDeltat = 0.5)
identical(round(gridDensity(theDens, pars, truncationPoint = 10, maxProbability = 0.01, maxDeltat = 0.5)$prob, digits = 10)
               ,round(gridDensity_Rcpp(theDens, pars, truncationPoint = 10, maxProbability = 0.01, maxDeltat = 0.5)$prob, digits = 10))

microbenchmark(gridDensity(theDens, pars, truncationPoint = 10, maxProbability = 0.01, maxDeltat = 0.5)
,gridDensity_Rcpp(theDens, pars, truncationPoint = 10, maxProbability = 0.01, maxDeltat = 0.5), times=1000L )

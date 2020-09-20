#include <Rcpp.h>
using namespace Rcpp;
#include "matrix_functions.h"

// Distributions

//' Phase-type density
//' 
//' Computes the density of phase-type distribution with parameters \code{pi} and \code{T} at \code{x}
//' @param x non-negative value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @return The density at \code{x}
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' phdensity(0.5, alpha, T) 
// [[Rcpp::export]]
NumericVector phdensity(NumericVector x, NumericVector pi, NumericMatrix T) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1), m_e);
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      density[k] = (1.0 - matrix_product(m_pi, m_e)(0,0));
    }
    else {
      density[k] = (matrix_product(m_pi, matrix_product(matrix_exponential(T * x[k]), m_t))(0,0));
    }
  }
  return density;
}


//' Phase-type cdf or tail
//' 
//' Computes the cdf of phase-type distribution with parameters \code{pi} and \code{T} at \code{x}
//' @param x non-negative value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @return The cdf (tail) at \code{x}
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' phcdf(0.5, alpha, T) 
//' phcdf(0.5, alpha, T, FALSE) 
// [[Rcpp::export]]
NumericVector phcdf(NumericVector x, NumericVector pi, NumericMatrix T, bool lower_tail = true) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      cdf[k] = (1.0 - matrix_product(m_pi, m_e)(0,0));
    }
    else {
      cdf[k] = (1.0 - matrix_product(m_pi, matrix_product(matrix_exponential(T * x[k]), m_e))(0,0));
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}


//' k moment of a phase-type
//' 
//' Computes the k moment of phase-type distribution with parameters \code{pi} and \code{T}
//' @param k Integer value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @return The k moment
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' phmoment(2, alpha, T) 
//' phmoment(4, alpha, T) 
// [[Rcpp::export]]
NumericVector phmoment(IntegerVector k, NumericVector pi, NumericMatrix T) {
  
  NumericVector moments(k.size());
  
  NumericVector fact_k = factorial(k);
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  
  for (int i = 0; i < k.size(); ++i){
    NumericMatrix U(matrix_inverse(T * (-1.0)));
    
    U = matrix_power(k[i], U);
    
    moments[i] = fact_k[i] * matrix_product(m_pi, matrix_product(U, m_e))(0,0);
  }
 return moments;
}

//' Laplace transform of a phase-type
//' 
//' Computes the Laplace transform at \code{s} of a phase-type distribution with parameters \code{pi} and \code{T}
//' @param s real value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @return Laplace transform
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' phLaplace(0.5, alpha, T) 
//' phLaplace(2.5, alpha, T) 
// [[Rcpp::export]]
NumericVector phLaplace(NumericVector s, NumericVector pi, NumericMatrix T) {
  
  NumericVector Laplace(s.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1), m_e);
  
  NumericMatrix identity_matrix = NumericMatrix::diag(T.nrow(), 1.0);
  
  for (int i = 0; i < s.size(); ++i){
    Laplace[i] = matrix_product(m_pi, matrix_product(matrix_inverse(matrix_sum(identity_matrix * s[i], T * (-1.0))), m_t))(0,0);
  }
  return Laplace;
}




//' IPH density - Slower
//' 
//' Computes the density of an IPH distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
//' @param x non-negative value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param g Tranformation 
//' @param g_inv Inverse of the transformation
//' @param lambda Derivative of the inverse
//' @param beta parameter of the transformation
//' @return The density at \code{x}
//' @examples
//' g <- function(x, beta) { x^(1/beta) }
//' g_inv <- function(x, beta) { x^beta}
//' lambda <- function(x, beta) {beta * x^(beta - 1)}
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' beta <- 0.5
//' iphdensity(0.5, alpha, T, g, g_inv, lambda, beta) 
// [[Rcpp::export]]
NumericVector iphdensity(NumericVector x, NumericVector pi, NumericMatrix T, Function g, Function g_inv, Function lambda, NumericVector beta) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1), m_e);
  
  NumericVector g_val = g(0, beta);
  NumericVector g_inv_val;
  NumericVector lambda_val;
  
  for (int k = 0; k < x.size(); ++k){
    g_inv_val = g_inv(x[k], beta);
    lambda_val = lambda(x[k], beta);
    if (x[k] == g_val[0]) {
      density[k] = (1.0 - matrix_product(m_pi, m_e)(0,0));
    }
    else {
      density[k] = (matrix_product(m_pi, matrix_product(matrix_exponential(T * g_inv_val[0]), m_t))(0,0)) * lambda_val[0];
    }
  }
  return density;
}



//' IPH cdf (tail)
//' 
//' Computes the cdf(tail) of an IPH distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
//' @param x non-negative value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param g Tranformation 
//' @param g_inv Inverse of the transformation
//' @param lambda Derivative of the inverse
//' @param beta parameter of the transformation
//' @return The cdf (tail) at \code{x}
//' @examples
//' g <- function(x, beta) { x^(1/beta) }
//' g_inv <- function(x, beta) { x^beta}
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' beta <- 0.5
//' iphcdf(0.5, alpha, T, g, g_inv, beta)
//' iphcdf(0.5, alpha, T, g, g_inv, beta, FALSE) 
// [[Rcpp::export]]
NumericVector iphcdf(NumericVector x, NumericVector pi, NumericMatrix T, Function g, Function g_inv, NumericVector beta, bool lower_tail = true) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  
  NumericVector g_val = g(0, beta);
  NumericVector g_inv_val;
  
  for (int k = 0; k < x.size(); ++k){
    g_inv_val = g_inv(x[k], beta);
    if (x[k] == g_val[0]) {
      cdf[k] = (1.0 - matrix_product(m_pi, m_e)(0,0));
    }
    else {
      cdf[k] = (1.0 - matrix_product(m_pi, matrix_product(matrix_exponential(T * g_inv_val[0]), m_e))(0,0));
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}



//' Matrix Weibull density
//' 
//' Computes the density of a matrix Weibull distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
//' @param x non-negative value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param beta shape parameter
//' @return The density at \code{x}
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' beta <- 0.5
//' mweibullden(0.5, alpha, T, beta) 
// [[Rcpp::export]]
NumericVector mWeibullden(NumericVector x, NumericVector pi, NumericMatrix T, double beta) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1), m_e);
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      density[k] = (1.0 - matrix_product(m_pi, m_e)(0,0));
    }
    else {
      density[k] = (matrix_product(m_pi, matrix_product(matrix_exponential(T * pow(x[k], beta)), m_t))(0,0)) * beta * pow(x[k], beta - 1);
    }
  }
  return density;
}



//' Matrix Weibull cdf
//' 
//' Computes the cdf (tail) of a matrix Weibull distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
//' @param x non-negative value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param beta shape parameter
//' @return The cdf (tail) at \code{x}
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' beta <- 0.5
//' mweibullcdf(0.5, alpha, T, beta) 
//' mweibullcdf(0.5, alpha, T, beta, FALSE) 
// [[Rcpp::export]]
NumericVector mWeibullcdf(NumericVector x, NumericVector pi, NumericMatrix T, double beta, bool lower_tail = true) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      cdf[k] = (1.0 - matrix_product(m_pi, m_e)(0,0));
    }
    else {
      cdf[k] = (1.0 - matrix_product(m_pi, matrix_product(matrix_exponential(T * pow(x[k], beta)), m_e))(0,0));
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}

// [[Rcpp::export]]
NumericVector RunFunction(NumericVector a, Function func)
{
  NumericVector b = func(a);
  return b;
}



//' Matrix Pareto density
//' 
//' Computes the density of a matrix Pareto distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
//' @param x non-negative value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param beta scale parameter
//' @return The density at \code{x}
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' beta <- 0.5
//' mparetoden(0.5, alpha, T, beta) 
// [[Rcpp::export]]
NumericVector mParetoden(NumericVector x, NumericVector pi, NumericMatrix T, double beta) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1), m_e);
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      density[k] = (1.0 - matrix_product(m_pi, m_e)(0,0));
    }
    else {
      density[k] = (matrix_product(m_pi, matrix_product(matrix_exponential(T * log(x[k] / beta + 1)), m_t))(0,0)) / (x[k] + beta);
    }
  }
  return density;
}



//' Matrix Pareto cdf
//' 
//' Computes the cdf (tail) of a matrix Pareto distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
//' @param x non-negative value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param beta shape parameter
//' @return The cdf (tail) at \code{x}
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' beta <- 0.5
//' mparetocdf(0.5, alpha, T, beta) 
//' mparetocdf(0.5, alpha, T, beta, FALSE) 
// [[Rcpp::export]]
NumericVector mParetocdf(NumericVector x, NumericVector pi, NumericMatrix T, double beta, bool lower_tail = true) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      cdf[k] = (1.0 - matrix_product(m_pi, m_e)(0,0));
    }
    else {
      cdf[k] = (1.0 - matrix_product(m_pi, matrix_product(matrix_exponential(T * log(x[k] / beta + 1)), m_e))(0,0));
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}



//' Matrix Gompertz density
//' 
//' Computes the density of a matrix Gompertz distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
//' @param x non-negative value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param beta  parameter
//' @return The density at \code{x}
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' beta <- 0.5
//' mgompertzden(0.5, alpha, T, beta) 
// [[Rcpp::export]]
NumericVector mGompertzden(NumericVector x, NumericVector pi, NumericMatrix T, double beta) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1), m_e);
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      density[k] = (1.0 - matrix_product(m_pi, m_e)(0,0));
    }
    else {
      density[k] = (matrix_product(m_pi, matrix_product(matrix_exponential(T * ((exp(x[k] * beta) - 1) / beta) ), m_t))(0,0)) * exp(x[k] * beta);
    }
  }
  return density;
}



//' Matrix Gompertz cdf
//' 
//' Computes the cdf (tail) of a matrix Gompertz distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
//' @param x non-negative value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param beta shape parameter
//' @return The cdf (tail) at \code{x}
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' beta <- 0.5
//' mgompertzcdf(0.5, alpha, T, beta) 
//' mgompertzcdf(0.5, alpha, T, beta, FALSE) 
// [[Rcpp::export]]
NumericVector mGompertzcdf(NumericVector x, NumericVector pi, NumericMatrix T, double beta, bool lower_tail = true) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      cdf[k] = (1.0 - matrix_product(m_pi, m_e)(0,0));
    }
    else {
      cdf[k] = (1.0 - matrix_product(m_pi, matrix_product(matrix_exponential(T * ((exp(x[k] * beta) - 1) / beta)), m_e))(0,0));
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}



//' Matrix GEV density
//' 
//' Computes the density of a matrix GEV distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
//' Dont allow for atoms in zero
//' @param x non-negative value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param beta  parameter
//' @return The density at \code{x}
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' mu <- 1
//' sigma <- 2
//' xi <- 0.5
//' mGEVden(0.5, alpha, T, mu, sigma, xi) 
// [[Rcpp::export]]
NumericVector mGEVDden(NumericVector x, NumericVector pi, NumericMatrix T, double mu, double sigma, double xi) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1), m_e);
  
  for (int k = 0; k < x.size(); ++k){
    if (xi == 0) {
      density[k] = (matrix_product(m_pi, matrix_product(matrix_exponential(T * exp(-(x[k] - mu) / sigma) ), m_t))(0,0)) * exp(-(x[k] - mu) / sigma) / sigma;
    }
    else {
      density[k] = (matrix_product(m_pi, matrix_product(matrix_exponential(T * pow(1 + (xi / sigma) * (x[k] - mu), -1 / xi) ), m_t))(0,0)) * pow(1 + (xi / sigma) * (x[k] - mu), -(1 + xi) / xi) / sigma;
    }
  }
  return density;
}



//' Matrix GEV cdf
//' 
//' Computes the cdf (tail) of a matrix GEV distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
//' @param x non-negative value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param beta shape parameter
//' @return The cdf (tail) at \code{x}
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' mu <- 1
//' sigma <- 2
//' xi <- 0.5
//' mGEVcdf(0.5, alpha, T, mu, sigma, xi) 
//' mGEVcdf(0.5, alpha, T, mu, sigma, xi, FALSE) 
// [[Rcpp::export]]
NumericVector mGEVDcdf(NumericVector x, NumericVector pi, NumericMatrix T, double mu, double sigma, double xi, bool lower_tail = true) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  
  for (int k = 0; k < x.size(); ++k){
    if (xi == 0) {
      cdf[k] = matrix_product(m_pi, matrix_product(matrix_exponential(T * exp(-(x[k] - mu) / sigma)), m_e))(0,0);
    }
    else {
      cdf[k] = matrix_product(m_pi, matrix_product(matrix_exponential(T * pow(1 + (xi / sigma) * (x[k] - mu), -1 / xi) ), m_e))(0,0);
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}


//' Bivariate phase-type joint density
//' 
//' @examples
//' alpha <- c(0.15, 0.85)
//' T11 <- matrix(c(c(-2,9),c(0,-11)), nrow = 2, ncol = 2)
//' T12 <- matrix(c(c(2,0),c(0,2)), nrow = 2, ncol = 2)
//' T22 <- matrix(c(c(-1,0),c(0.5,-5)), nrow = 2, ncol = 2)
//' x1 <- matrix(c(0.5,2), ncol=2) 
//' x2 <- matrix(c(c(0.5,1), c(2, 1.5)), ncol=2) 
//' bivphden(x1, alpha, T11, T12, T22) 
//' bivphden(x2, alpha, T11, T12, T22) 
// [[Rcpp::export]]
NumericVector bivphden(NumericMatrix x, NumericVector alpha, NumericMatrix T11, NumericMatrix T12, NumericMatrix T22) {
  
  long N{x.nrow()};
  long p2{T22.nrow()};
  
  NumericVector density(N);
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(p2, 1);
  NumericMatrix m_e(p2, 1, e.begin());
  NumericMatrix m_t = matrix_product(T22 * (-1), m_e);
  
  for (int k = 0; k < N; ++k){
   density[k] = matrix_product(m_alpha, matrix_product(matrix_exponential(T11 * x(k,0)), matrix_product(T12, matrix_product(matrix_exponential(T22 * x(k,1)), m_t))))(0,0);
  }
  
  return density;
}


//' Bivariate phase-type joint tail
//' 
//' @examples
//' alpha <- c(0.15, 0.85)
//' T11 <- matrix(c(c(-2,9),c(0,-11)), nrow = 2, ncol = 2)
//' T12 <- matrix(c(c(2,0),c(0,2)), nrow = 2, ncol = 2)
//' T22 <- matrix(c(c(-1,0),c(0.5,-5)), nrow = 2, ncol = 2)
//' x1 <- matrix(c(0.5,1), ncol=2) 
//' x2 <- matrix(c(c(0.5,1), c(2, 1.5)), ncol=2) 
//' bivphtail(x1, alpha, T11, T12, T22) 
//' bivphtail(x2, alpha, T11, T12, T22) 
// [[Rcpp::export]]
NumericVector bivphtail(NumericMatrix x, NumericVector alpha, NumericMatrix T11, NumericMatrix T12, NumericMatrix T22) {
  
  long N{x.nrow()};
  long p2{T22.nrow()};
  
  NumericVector tail(N);
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(p2, 1);
  NumericMatrix m_e(p2, 1, e.begin());
  NumericMatrix m_t = matrix_product(T22 * (-1), m_e);
  
  for (int k = 0; k < N; ++k){
    tail[k] = matrix_product(m_alpha, matrix_product(matrix_inverse(T11 * (-1.0)), matrix_product(matrix_exponential(T11 * x(k,0)), matrix_product(T12, matrix_product( matrix_exponential(T22 * x(k,1)), m_e)))))(0,0);
  }
  
  return tail;
}



//' Bivariate matrix Weibull joint density
//' 
//' @examples
//' alpha <- c(0.15, 0.85)
//' T11 <- matrix(c(c(-2,9),c(0,-11)), nrow = 2, ncol = 2)
//' T12 <- matrix(c(c(2,0),c(0,2)), nrow = 2, ncol = 2)
//' T22 <- matrix(c(c(-1,0),c(0.5,-5)), nrow = 2, ncol = 2)
//' beta <- c(0.5, 0.7)
//' x1 <- matrix(c(0.5,2), ncol=2) 
//' x2 <- matrix(c(c(0.5,1), c(2, 1.5)), ncol=2) 
//' bivmWeibden(x1, alpha, T11, T12, T22, beta) 
//' bivmWeibden(x2, alpha, T11, T12, T22, beta) 
// [[Rcpp::export]]
NumericVector bivmWeibden(NumericMatrix x, NumericVector alpha, NumericMatrix T11, NumericMatrix T12, NumericMatrix T22, NumericVector beta) {
  
  long N{x.nrow()};
  long p2{T22.nrow()};
  
  NumericVector density(N);
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(p2, 1);
  NumericMatrix m_e(p2, 1, e.begin());
  NumericMatrix m_t = matrix_product(T22 * (-1), m_e);
  
  for (int k = 0; k < N; ++k){
    density[k] = matrix_product(m_alpha, matrix_product(matrix_exponential(T11 * pow(x(k,0), beta[0])), matrix_product(T12, matrix_product(matrix_exponential(T22 * pow(x(k,1), beta[1])), m_t))))(0,0) * beta[0] *  beta[1] * pow(x(k,0), beta[0] - 1) * pow(x(k,1), beta[1] - 1);
  }
  
  return density;
}


//' Bivariate matrix Weibull joint tail
//' 
//' @examples
//' alpha <- c(0.15, 0.85)
//' T11 <- matrix(c(c(-2,9),c(0,-11)), nrow = 2, ncol = 2)
//' T12 <- matrix(c(c(2,0),c(0,2)), nrow = 2, ncol = 2)
//' T22 <- matrix(c(c(-1,0),c(0.5,-5)), nrow = 2, ncol = 2)
//' beta <- c(0.5, 0.7)
//' x1 <- matrix(c(0.5,1), ncol=2) 
//' x2 <- matrix(c(c(0.5,1), c(2, 1.5)), ncol=2) 
//' bimWeibtail(x1, alpha, T11, T12, T22, beta) 
//' bimWeibtail(x2, alpha, T11, T12, T22, beta) 
// [[Rcpp::export]]
NumericVector bimWeibtail(NumericMatrix x, NumericVector alpha, NumericMatrix T11, NumericMatrix T12, NumericMatrix T22, NumericVector beta) {
  
  long N{x.nrow()};
  long p2{T22.nrow()};
  
  NumericVector tail(N);
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(p2, 1);
  NumericMatrix m_e(p2, 1, e.begin());
  NumericMatrix m_t = matrix_product(T22 * (-1), m_e);
  
  for (int k = 0; k < N; ++k){
    tail[k] = matrix_product(m_alpha, matrix_product(matrix_inverse(T11 * (-1.0)), matrix_product(matrix_exponential(T11 * pow(x(k,0), beta[0])), matrix_product(T12, matrix_product( matrix_exponential(T22 * pow(x(k,1), beta[1])), m_e)))))(0,0);
  }
  
  return tail;
}


//' Bivariate matrix Pareto joint density
//' 
//' @examples
//' alpha <- c(0.15, 0.85)
//' T11 <- matrix(c(c(-2,9),c(0,-11)), nrow = 2, ncol = 2)
//' T12 <- matrix(c(c(2,0),c(0,2)), nrow = 2, ncol = 2)
//' T22 <- matrix(c(c(-1,0),c(0.5,-5)), nrow = 2, ncol = 2)
//' beta <- c(2, 4)
//' x1 <- matrix(c(0.5,2), ncol=2) 
//' x2 <- matrix(c(c(0.5,1), c(2, 1.5)), ncol=2) 
//' bivmParden(x1, alpha, T11, T12, T22, beta) 
//' bivmParden(x2, alpha, T11, T12, T22, beta) 
// [[Rcpp::export]]
NumericVector bivmParden(NumericMatrix x, NumericVector alpha, NumericMatrix T11, NumericMatrix T12, NumericMatrix T22, NumericVector beta) {
  
  long N{x.nrow()};
  long p2{T22.nrow()};
  
  NumericVector density(N);
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(p2, 1);
  NumericMatrix m_e(p2, 1, e.begin());
  NumericMatrix m_t = matrix_product(T22 * (-1), m_e);
  
  for (int k = 0; k < N; ++k){
    density[k] = matrix_product(m_alpha, matrix_product(matrix_exponential(T11 * log(x(k,0) / beta[0] + 1)), matrix_product(T12, matrix_product(matrix_exponential(T22 * log(x(k,1) / beta[1] + 1)), m_t))))(0,0) / ((x(k,0) + beta[0]) * (x(k,1) + beta[1]));
  }
  
  return density;
}



//' Bivariate matrix Weibull joint tail
//' 
//' @examples
//' alpha <- c(0.15, 0.85)
//' T11 <- matrix(c(c(-2,9),c(0,-11)), nrow = 2, ncol = 2)
//' T12 <- matrix(c(c(2,0),c(0,2)), nrow = 2, ncol = 2)
//' T22 <- matrix(c(c(-1,0),c(0.5,-5)), nrow = 2, ncol = 2)
//' beta <- c(2, 4)
//' x1 <- matrix(c(0.5,1), ncol=2) 
//' x2 <- matrix(c(c(0.5,1), c(2, 1.5)), ncol=2) 
//' bimPartail(x1, alpha, T11, T12, T22, beta) 
//' bimPartail(x2, alpha, T11, T12, T22, beta) 
// [[Rcpp::export]]
NumericVector bimPartail(NumericMatrix x, NumericVector alpha, NumericMatrix T11, NumericMatrix T12, NumericMatrix T22, NumericVector beta) {
  
  long N{x.nrow()};
  long p2{T22.nrow()};
  
  NumericVector tail(N);
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(p2, 1);
  NumericMatrix m_e(p2, 1, e.begin());
  NumericMatrix m_t = matrix_product(T22 * (-1), m_e);
  
  for (int k = 0; k < N; ++k){
    tail[k] = matrix_product(m_alpha, matrix_product(matrix_inverse(T11 * (-1.0)), matrix_product(matrix_exponential(T11 * log(x(k,0) / beta[0] + 1)), matrix_product(T12, matrix_product( matrix_exponential(T22 * log(x(k,1) / beta[1] + 1)), m_e)))))(0,0);
  }
  
  return tail;
}


//' Pi and T of a linear combination of a MPH*
//' 
//' @examples
//' pi <- c(0.15, 0.85, 0 ,0)
//' T11 <- matrix(c(c(-2,9),c(0,-11)), nrow = 2, ncol = 2)
//' T12 <- matrix(c(c(2,0),c(0,2)), nrow = 2, ncol = 2)
//' T22 <- matrix(c(c(-1,0),c(0.5,-5)), nrow = 2, ncol = 2)
//' T <- merge_matrices(T11, T12, T22)
//' R <- matrix(c(c(1,1,0,0), c(0,0,1,1)), ncol=2)
//' w1 <- c(1,0)
//' linear_combination(w1, pi, T, R)
//' w2 <- c(0,1)
//' linear_combination(w2, pi, T, R)
//' matrix(c(0.15, 0.85), ncol=2)%*%matrix_inverse(T11 * (-1))%*%T12
//' w3 <- c(1,1)
//' linear_combination(w3, pi, T, R)
// [[Rcpp::export]]
List linear_combination(NumericVector w, NumericVector pi, NumericMatrix T, NumericMatrix R) {
  long p{T.nrow()};
  
  NumericVector newStates;
  
  int NumZeros{0};
  IntegerVector deleteRows; //states to be deleted
  IntegerVector keepRows; //states to keep
  
  NumericMatrix m_w(w.size(), 1, w.begin());
  
  NumericMatrix Rw(matrix_product(R, m_w)); 
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  
  for (int j{0}; j < p; ++j) {
    if (Rw(j,0) == 0) {
      deleteRows.push_back(j);
      ++NumZeros;
    }
    else {
      keepRows.push_back(j);
    }
  }
  
  newStates = keepRows;
  NumericVector pi_w(p - NumZeros);
  NumericMatrix T_w(p - NumZeros, p - NumZeros);
  
  if (NumZeros == 0) {
    NumericMatrix diagonal(p,p);
    
    for (int i{0}; i < p; ++i) {
      diagonal(i,i) = 1.0 / Rw(i,0);
    }
    diagonal = matrix_product(diagonal, T);
    
    pi_w = pi;
    T_w = diagonal;
    
  }
  else {
    long n1{deleteRows.size()};
    long n2{keepRows.size()};
    
    NumericMatrix Spp(n2,n2);
    NumericMatrix Sp0(n2,n1);
    NumericMatrix S0p(n1,n2);
    NumericMatrix S00(n1,n1);
    
    NumericMatrix Taux(n2,n2);
    NumericMatrix diagonal(n2,n2);
    
    NumericMatrix pi0(1,n1);
    NumericMatrix pip(1,n2);
    
    NumericMatrix piaux(1,n2);
    
    for (int i{0}; i < n2; i++) {
      for (int j = 0; j < n2; j++) {
        Spp(i,j) = T(keepRows[i],keepRows[j]);
      }
      for (int j{0}; j < n1; j++) {
        Sp0(i,j) = T(keepRows[i],deleteRows[j]);
      }
      pip(0,i) = m_pi(0,keepRows[i]);
    }
    for (int i{0}; i < n1; i++) {
      for (int j{0}; j < n2; j++) {
        S0p(i,j) = T(deleteRows[i],keepRows[j]);
      }
      for (int j{0}; j < n1; j++){
        S00(i,j) = T(deleteRows[i],deleteRows[j]);
      }
      pi0(0,i) = m_pi(0,deleteRows[i]);
    }
    
    piaux = matrix_sum(pip, matrix_product(pi0, matrix_product(matrix_inverse(S00 * (-1.0)), S0p)));
    
    Taux = matrix_sum(Spp, matrix_product(Sp0, matrix_product(matrix_inverse(S00 * (-1.0)), S0p)));
    
    for (int i{0}; i < n2; ++i) {
      diagonal(i,i) = 1.0 / Rw(keepRows[i],0);
      pi_w[i] = piaux(0, i);
    }
    diagonal = matrix_product(diagonal, Taux);
    
    T_w = diagonal;
    
  }
  keepRows = keepRows + 1;
  List L = List::create(Named("piw") = pi_w, _["Tw"] = T_w, _["new_states"] = keepRows);
  
  return L;
}

//' Joint MGF of a MPH
// [[Rcpp::export]]
double jointMGF(const NumericVector & w, NumericVector pi, NumericMatrix T, NumericMatrix R) {
  long p{T.nrow()};
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(p, 1);
  NumericMatrix m_e(p, 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1), m_e);
  
  NumericMatrix m_w(w.size(), 1, w.begin());
  
  NumericMatrix Rw(matrix_product(R, m_w));
  
  NumericMatrix diagonal(p,p);
  
  for (int i{0}; i < p; ++i) {
    diagonal(i,i) = -Rw(i,0);
  }
  return matrix_product(m_pi, matrix_product(matrix_inverse(matrix_sum(diagonal, T * (-1.0))), m_t))(0,0);
}



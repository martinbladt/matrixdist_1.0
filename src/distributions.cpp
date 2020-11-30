#include <Rcpp.h>
using namespace Rcpp;
#include "matrix_functions.h"
#include "distributions.h"

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
NumericVector mWeibullcdf(NumericVector x, NumericVector pi, NumericMatrix T, double beta, bool lower_tail) {
  
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


//' Matrix LogNormal density
//' 
//' Computes the density of a matrix LogNormal distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
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
NumericVector mLogNormalden(NumericVector x, NumericVector pi, NumericMatrix T, double beta) {
  
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
      density[k] = (matrix_product(m_pi, matrix_product(matrix_exponential(T * pow(log(x[k] + 1), beta)), m_t))(0,0)) * beta * pow(log(x[k] + 1), beta - 1)/(x[k] + 1);
    }
  }
  return density;
}



//' Matrix LogNormal cdf
//' 
//' Computes the cdf (tail) of a matrix LogNormal distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
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
NumericVector mLogNormalcdf(NumericVector x, NumericVector pi, NumericMatrix T, double beta, bool lower_tail = true) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      cdf[k] = (1.0 - matrix_product(m_pi, m_e)(0,0));
    }
    else {
      cdf[k] = (1.0 - matrix_product(m_pi, matrix_product(matrix_exponential(T * pow(log(x[k] + 1), beta)), m_e))(0,0));
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}



//' Matrix Log-Logistic density
//' 
//' Computes the density of a matrix Log-Logistic distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
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
NumericVector mLogLogisticden(NumericVector x, NumericVector pi, NumericMatrix T, NumericVector beta) {
  
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
      density[k] = (matrix_product(m_pi, matrix_product(matrix_exponential(T * log(pow(x[k] / beta[0], beta[1]) + 1)), m_t))(0,0)) * (pow(x[k] / beta[0], beta[1] - 1) * beta[1] / beta[0]) / (pow(x[k] / beta[0], beta[1]) + 1);
    }
  }
  return density;
}


//' Matrix Log-Logistic cdf
//' 
//' Computes the cdf (tail) of a matrix Log-Logistic distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
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
NumericVector mLogLogisticcdf(NumericVector x, NumericVector pi, NumericMatrix T, NumericVector beta, bool lower_tail = true) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      cdf[k] = (1.0 - matrix_product(m_pi, m_e)(0,0));
    }
    else {
      cdf[k] = (1.0 - matrix_product(m_pi, matrix_product(matrix_exponential(T * log(pow(x[k] / beta[0], beta[1]) + 1)), m_e))(0,0));
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
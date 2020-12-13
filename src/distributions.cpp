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
//' 
// [[Rcpp::export]]
NumericVector phdensity(NumericVector x, NumericVector pi, NumericMatrix T) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1.0), m_e);
  
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
//' @param lower_tail cdf or tail
//' @return The cdf (tail) at \code{x}
//' 
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
//' 
// [[Rcpp::export]]
NumericVector phmoment(NumericVector k, NumericVector pi, NumericMatrix T) {
  
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
//' 
// [[Rcpp::export]]
NumericVector phLaplace(NumericVector s, NumericVector pi, NumericMatrix T) {
  
  NumericVector Laplace(s.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1.0), m_e);
  
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
//' 
// [[Rcpp::export]]
NumericVector mweibullden(NumericVector x, NumericVector pi, NumericMatrix T, double beta) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1.0), m_e);
  
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
//' @param lower_tail cdf or tail
//' @return The cdf (tail) at \code{x}
//' 
// [[Rcpp::export]]
NumericVector mweibullcdf(NumericVector x, NumericVector pi, NumericMatrix T, double beta, bool lower_tail) {
  
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
    return (1.0 - cdf);
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
//' 
// [[Rcpp::export]]
NumericVector mparetoden(NumericVector x, NumericVector pi, NumericMatrix T, double beta) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1.0), m_e);
  
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
//' @param lower_tail cdf or tail
//' @return The cdf (tail) at \code{x}
//' 
// [[Rcpp::export]]
NumericVector mparetocdf(NumericVector x, NumericVector pi, NumericMatrix T, double beta, bool lower_tail = true) {
  
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
    return (1.0 - cdf);
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
//' 
// [[Rcpp::export]]
NumericVector mlognormalden(NumericVector x, NumericVector pi, NumericMatrix T, double beta) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1.0), m_e);
  
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
//' @param lower_tail cdf or tail
//' @return The cdf (tail) at \code{x}
//' 
// [[Rcpp::export]]
NumericVector mlognormalcdf(NumericVector x, NumericVector pi, NumericMatrix T, double beta, bool lower_tail = true) {
  
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
    return (1.0 - cdf);
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
//' 
// [[Rcpp::export]]
NumericVector mloglogisticden(NumericVector x, NumericVector pi, NumericMatrix T, NumericVector beta) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1.0), m_e);
  
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
//' @param lower_tail cdf or tail
//' @return The cdf (tail) at \code{x}
//' 
// [[Rcpp::export]]
NumericVector mloglogisticcdf(NumericVector x, NumericVector pi, NumericMatrix T, NumericVector beta, bool lower_tail = true) {
  
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
    return (1.0 - cdf);
  }
}


//' Matrix Gompertz density
//' 
//' Computes the density of a matrix Gompertz distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
//' @param x non-negative value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param beta parameter
//' @return The density at \code{x}
//' 
// [[Rcpp::export]]
NumericVector mgompertzden(NumericVector x, NumericVector pi, NumericMatrix T, double beta) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1.0), m_e);
  
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
//' @param lower_tail cdf or tail
//' @return The cdf (tail) at \code{x}
//' 
// [[Rcpp::export]]
NumericVector mgompertzcdf(NumericVector x, NumericVector pi, NumericMatrix T, double beta, bool lower_tail = true) {
  
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
    return (1.0 - cdf);
  }
}



//' Matrix GEV density
//' 
//' Computes the density of a matrix GEV distribution with parameters \code{pi}, \code{T} and \code{beta} at \code{x}
//' Dont allow for atoms in zero
//' @param x non-negative value
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param mu  location parameter
//' @param sigma scale parameter
//' @param xi shape parameter
//' @return The density at \code{x}
//' 
// [[Rcpp::export]]
NumericVector mgevden(NumericVector x, NumericVector pi, NumericMatrix T, double mu, double sigma, double xi) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1.0), m_e);
  
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
//' @param mu  location parameter
//' @param sigma scale parameter
//' @param xi shape parameter
//' @param lower_tail cdf or tail
//' @return The cdf (tail) at \code{x}
//' 
// [[Rcpp::export]]
NumericVector mgevcdf(NumericVector x, NumericVector pi, NumericMatrix T, double mu, double sigma, double xi, bool lower_tail = true) {
  
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
    return (1.0 - cdf);
  }
}
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




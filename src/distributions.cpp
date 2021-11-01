#include <Rcpp.h>
using namespace Rcpp;
#include "exp_arm.h"

// Distributions

//' Phase-type density
//' 
//' Computes the density of phase-type distribution with parameters \code{alpha}
//'  and \code{S} at \code{x}.
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
NumericVector phdensity(NumericVector x, NumericVector alpha, NumericMatrix S) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(alpha.size(), 1.0);
  NumericMatrix m_e(alpha.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(S * (-1.0), m_e);
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      density[k] = (1.0 - matrix_product(m_alpha, m_e)(0,0));
    }
    else {
      density[k] = (matrix_product(m_alpha, matrix_product(matrix_exponential(S * x[k]), m_t))(0,0));
    }
  }
  return density;
}


//' Phase-type cdf or tail
//' 
//' Computes the cdf of phase-type distribution with parameters \code{alpha} and
//'  \code{S} at \code{x}.
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}
//' 
// [[Rcpp::export]]
NumericVector phcdf(NumericVector x, NumericVector alpha, NumericMatrix S, bool lower_tail = true) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(alpha.size(), 1.0);
  NumericMatrix m_e(alpha.size(), 1, e.begin());
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      cdf[k] = (1.0 - matrix_product(m_alpha, m_e)(0,0));
    }
    else {
      cdf[k] = (1.0 - matrix_product(m_alpha, matrix_product(matrix_exponential(S * x[k]), m_e))(0,0));
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
//' Computes the density of a matrix Weibull distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Shape parameter.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
NumericVector mweibullden(NumericVector x, NumericVector alpha, NumericMatrix S, double beta) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(alpha.size(), 1.0);
  NumericMatrix m_e(alpha.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(S * (-1.0), m_e);
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      density[k] = (1.0 - matrix_product(m_alpha, m_e)(0,0));
    }
    else {
      density[k] = (matrix_product(m_alpha, matrix_product(matrix_exponential(S * pow(x[k], beta)), m_t))(0,0)) * beta * pow(x[k], beta - 1.0);
    }
  }
  return density;
}



//' Matrix Weibull cdf
//' 
//' Computes the cdf (tail) of a matrix Weibull distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Shape parameter.
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}.
//' 
// [[Rcpp::export]]
NumericVector mweibullcdf(NumericVector x, NumericVector alpha, NumericMatrix S, double beta, bool lower_tail) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(alpha.size(), 1.0);
  NumericMatrix m_e(alpha.size(), 1, e.begin());
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      cdf[k] = (1.0 - matrix_product(m_alpha, m_e)(0,0));
    }
    else {
      cdf[k] = (1.0 - matrix_product(m_alpha, matrix_product(matrix_exponential(S * pow(x[k], beta)), m_e))(0,0));
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1.0 - cdf);
  }
}


//' Matrix Pareto density
//' 
//' Computes the density of a matrix Pareto distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Scale parameter.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
NumericVector mparetoden(NumericVector x, NumericVector alpha, NumericMatrix S, double beta) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(alpha.size(), 1.0);
  NumericMatrix m_e(alpha.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(S * (-1.0), m_e);
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      density[k] = (1.0 - matrix_product(m_alpha, m_e)(0,0));
    }
    else {
      density[k] = (matrix_product(m_alpha, matrix_product(matrix_exponential(S * log(x[k] / beta + 1.0)), m_t))(0,0)) / (x[k] + beta);
    }
  }
  return density;
}



//' Matrix Pareto cdf
//' 
//' Computes the cdf (tail) of a matrix Pareto distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Shape parameter.
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}.
//' 
// [[Rcpp::export]]
NumericVector mparetocdf(NumericVector x, NumericVector alpha, NumericMatrix S, double beta, bool lower_tail = true) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(alpha.size(), 1.0);
  NumericMatrix m_e(alpha.size(), 1, e.begin());
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      cdf[k] = (1.0 - matrix_product(m_alpha, m_e)(0,0));
    }
    else {
      cdf[k] = (1.0 - matrix_product(m_alpha, matrix_product(matrix_exponential(S * log(x[k] / beta + 1.0)), m_e))(0,0));
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
//' Computes the density of a matrix LogNormal distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Shape parameter.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
NumericVector mlognormalden(NumericVector x, NumericVector alpha, NumericMatrix S, double beta) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(alpha.size(), 1.0);
  NumericMatrix m_e(alpha.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(S * (-1.0), m_e);
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      density[k] = (1.0 - matrix_product(m_alpha, m_e)(0,0));
    }
    else {
      density[k] = (matrix_product(m_alpha, matrix_product(matrix_exponential(S * pow(log(x[k] + 1.0), beta)), m_t))(0,0)) * beta * pow(log(x[k] + 1), beta - 1)/(x[k] + 1);
    }
  }
  return density;
}



//' Matrix LogNormal cdf
//' 
//' Computes the cdf (tail) of a matrix LogNormal distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Shape parameter.
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}.
//' 
// [[Rcpp::export]]
NumericVector mlognormalcdf(NumericVector x, NumericVector alpha, NumericMatrix S, double beta, bool lower_tail = true) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(alpha.size(), 1.0);
  NumericMatrix m_e(alpha.size(), 1, e.begin());
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      cdf[k] = (1.0 - matrix_product(m_alpha, m_e)(0,0));
    }
    else {
      cdf[k] = (1.0 - matrix_product(m_alpha, matrix_product(matrix_exponential(S * pow(log(x[k] + 1.0), beta)), m_e))(0,0));
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
//' Computes the density of a matrix Log-Logistic distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Scale parameter.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
NumericVector mloglogisticden(NumericVector x, NumericVector alpha, NumericMatrix S, NumericVector beta) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(alpha.size(), 1.0);
  NumericMatrix m_e(alpha.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(S * (-1.0), m_e);
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      density[k] = (1.0 - matrix_product(m_alpha, m_e)(0,0));
    }
    else {
      density[k] = (matrix_product(m_alpha, matrix_product(matrix_exponential(S * log(pow(x[k] / beta[0], beta[1]) + 1.0)), m_t))(0,0)) * (pow(x[k] / beta[0], beta[1] - 1) * beta[1] / beta[0]) / (pow(x[k] / beta[0], beta[1]) + 1);
    }
  }
  return density;
}


//' Matrix Log-Logistic cdf
//' 
//' Computes the cdf (tail) of a matrix Log-Logistic distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Shape parameter.
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}.
//' 
// [[Rcpp::export]]
NumericVector mloglogisticcdf(NumericVector x, NumericVector alpha, NumericMatrix S, NumericVector beta, bool lower_tail = true) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(alpha.size(), 1.0);
  NumericMatrix m_e(alpha.size(), 1, e.begin());
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      cdf[k] = (1.0 - matrix_product(m_alpha, m_e)(0,0));
    }
    else {
      cdf[k] = (1.0 - matrix_product(m_alpha, matrix_product(matrix_exponential(S * log(pow(x[k] / beta[0], beta[1]) + 1.0)), m_e))(0,0));
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
//' Computes the density of a matrix Gompertz distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Parameter.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
NumericVector mgompertzden(NumericVector x, NumericVector alpha, NumericMatrix S, double beta) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(alpha.size(), 1.0);
  NumericMatrix m_e(alpha.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(S * (-1.0), m_e);
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      density[k] = (1.0 - matrix_product(m_alpha, m_e)(0,0));
    }
    else {
      density[k] = (matrix_product(m_alpha, matrix_product(matrix_exponential(S * ((exp(x[k] * beta) - 1.0) / beta) ), m_t))(0,0)) * exp(x[k] * beta);
    }
  }
  return density;
}



//' Matrix Gompertz cdf
//' 
//' Computes the cdf (tail) of a matrix Gompertz distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Shape parameter.
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}.
//' 
// [[Rcpp::export]]
NumericVector mgompertzcdf(NumericVector x, NumericVector alpha, NumericMatrix S, double beta, bool lower_tail = true) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(alpha.size(), 1.0);
  NumericMatrix m_e(alpha.size(), 1, e.begin());
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      cdf[k] = (1.0 - matrix_product(m_alpha, m_e)(0,0));
    }
    else {
      cdf[k] = (1.0 - matrix_product(m_alpha, matrix_product(matrix_exponential(S * ((exp(x[k] * beta) - 1.0) / beta)), m_e))(0,0));
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
//' Computes the density of a matrix GEV distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' Dont allow for atoms in zero
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param mu Location parameter.
//' @param sigma Scale parameter.
//' @param xi Shape parameter.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
NumericVector mgevden(NumericVector x, NumericVector alpha, NumericMatrix S, double mu, double sigma, double xi) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(alpha.size(), 1.0);
  NumericMatrix m_e(alpha.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(S * (-1.0), m_e);
  
  for (int k = 0; k < x.size(); ++k){
    if (xi == 0) {
      density[k] = (matrix_product(m_alpha, matrix_product(matrix_exponential(S * exp(-(x[k] - mu) / sigma) ), m_t))(0,0)) * exp(-(x[k] - mu) / sigma) / sigma;
    }
    else {
      density[k] = (matrix_product(m_alpha, matrix_product(matrix_exponential(S * pow(1.0 + (xi / sigma) * (x[k] - mu), -1.0 / xi) ), m_t))(0,0)) * pow(1 + (xi / sigma) * (x[k] - mu), -(1 + xi) / xi) / sigma;
    }
  }
  return density;
}



//' Matrix GEV cdf
//' 
//' Computes the cdf (tail) of a matrix GEV distribution with parameters 
//' \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param mu Location parameter.
//' @param sigma Scale parameter.
//' @param xi Shape parameter.
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}.
//' 
// [[Rcpp::export]]
NumericVector mgevcdf(NumericVector x, NumericVector alpha, NumericMatrix S, double mu, double sigma, double xi, bool lower_tail = true) {
  
  NumericVector cdf(x.size());
  
  NumericMatrix m_alpha(1, alpha.size(), alpha.begin());
  NumericVector e(alpha.size(), 1.0);
  NumericMatrix m_e(alpha.size(), 1, e.begin());
  
  for (int k = 0; k < x.size(); ++k){
    if (xi == 0) {
      cdf[k] = matrix_product(m_alpha, matrix_product(matrix_exponential(S * exp(-(x[k] - mu) / sigma)), m_e))(0,0);
    }
    else {
      cdf[k] = matrix_product(m_alpha, matrix_product(matrix_exponential(S * pow(1.0 + (xi / sigma) * (x[k] - mu), -1.0 / xi) ), m_e))(0,0);
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1.0 - cdf);
  }
}
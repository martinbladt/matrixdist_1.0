# include <RcppArmadillo.h>
#include "m_exp.h"
// [[ Rcpp :: depends ( RcppArmadillo )]]

// Distributions

//' Phase-type density
//' 
//' Computes the density of a phase-type distribution with parameters
//'  \code{alpha} and \code{S} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector phdensity(Rcpp::NumericVector x, arma::vec alpha, arma::mat S) {
  Rcpp::NumericVector density(x.size());
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    if (x[k] == 0) {
      aux_mat = alpha.t() * e;
      density[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * x[k]) * exit_vect;
      density[k] = aux_mat(0,0);
    }
  }
  return density;
}


//' Phase-type cdf
//' 
//' Computes the cdf (tail) of a phase-type distribution with parameters 
//'  \code{alpha} and \code{S} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector phcdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, bool lower_tail = true) {
  Rcpp::NumericVector cdf(x.size());
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    if (x[k] == 0) {
      aux_mat = alpha.t() * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * x[k]) * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}


//' Matrix-Weibull density
//' 
//' Computes the density of a matrix-Weibull distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Shape parameter.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector mweibullden(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta) {
  Rcpp::NumericVector density(x.size());
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    if (x[k] == 0) {
      aux_mat = alpha.t() * e;
      density[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * pow(x[k], beta)) * exit_vect;
      density[k] = aux_mat(0,0) * beta * pow(x[k], beta - 1.0);
    }
  }
  return density;
}


//' Matrix-Weibull cdf
//' 
//' Computes the cdf (tail) of a matrix-Weibull distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Shape parameter.
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector mweibullcdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta, bool lower_tail = true) {
  Rcpp::NumericVector cdf(x.size());
  
  arma::mat e; 
  e.ones(S.n_cols, 1);
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    if (x[k] == 0) {
      aux_mat = alpha.t() * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * pow(x[k], beta)) * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}


//' Matrix-Pareto density
//' 
//' Computes the density of a matrix-Pareto distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Scale parameter.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector mparetoden(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta) {
  Rcpp::NumericVector density(x.size());
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    if (x[k] == 0) {
      aux_mat = alpha.t() * e;
      density[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * log(x[k] / beta + 1.0)) * exit_vect;
      density[k] = aux_mat(0,0) / (x[k] + beta);
    }
  }
  return density;
}


//' Matrix-Pareto cdf
//' 
//' Computes the cdf (tail) of a matrix-Pareto distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Scale parameter.
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector mparetocdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta, bool lower_tail = true) {
  Rcpp::NumericVector cdf(x.size());
  
  arma::mat e; 
  e.ones(S.n_cols, 1);
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    if (x[k] == 0) {
      aux_mat = alpha.t() * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * log(x[k] / beta + 1.0)) * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}


//' Matrix-lognormal density
//' 
//' Computes the density of a matrix-lognormal distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Shape parameter.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector mlognormalden(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta) {
  Rcpp::NumericVector density(x.size());
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    if (x[k] == 0) {
      aux_mat = alpha.t() * e;
      density[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * pow(log(x[k] + 1.0), beta)) * exit_vect;
      density[k] = aux_mat(0,0) * beta * pow(log(x[k] + 1), beta - 1)/(x[k] + 1);
    }
  }
  return density;
}


//' Matrix-lognormal cdf
//' 
//' Computes the cdf (tail) of a matrix-lognormal distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Shape parameter.
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector mlognormalcdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta, bool lower_tail = true) {
  Rcpp::NumericVector cdf(x.size());
  
  arma::mat e; 
  e.ones(S.n_cols, 1);
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    if (x[k] == 0) {
      aux_mat = alpha.t() * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * pow(log(x[k] + 1.0), beta)) * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}


//' Matrix-loglogistic density
//' 
//' Computes the density of a matrix-loglogistic distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Transformation parameters.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector mloglogisticden(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, Rcpp::NumericVector beta) {
  Rcpp::NumericVector density(x.size());
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    if (x[k] == 0) {
      aux_mat = alpha.t() * e;
      density[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * log(pow(x[k] / beta[0], beta[1]) + 1.0)) * exit_vect;
      density[k] = aux_mat(0,0) * (pow(x[k] / beta[0], beta[1] - 1) * beta[1] / beta[0]) / (pow(x[k] / beta[0], beta[1]) + 1);
    }
  }
  return density;
}


//' Matrix-loglogistic cdf
//' 
//' Computes the cdf (tail) of a matrix-loglogistic distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Transformation parameters.
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector mloglogisticcdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, Rcpp::NumericVector beta, bool lower_tail = true) {
  Rcpp::NumericVector cdf(x.size());
  
  arma::mat e; 
  e.ones(S.n_cols, 1);
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    if (x[k] == 0) {
      aux_mat = alpha.t() * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * log(pow(x[k] / beta[0], beta[1]) + 1.0)) * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}


//' Matrix-Gompertz density
//' 
//' Computes the density of a matrix-Gompertz distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Shape parameter.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector mgompertzden(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta) {
  Rcpp::NumericVector density(x.size());
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    if (x[k] == 0) {
      aux_mat = alpha.t() * e;
      density[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * ((exp(x[k] * beta) - 1.0) / beta)) * exit_vect;
      density[k] = aux_mat(0,0) * exp(x[k] * beta);
    }
  }
  return density;
}


//' Matrix-Gompertz cdf
//' 
//' Computes the cdf (tail) of a matrix-Gompertz distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Shape parameter.
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector mgompertzcdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta, bool lower_tail = true) {
  Rcpp::NumericVector cdf(x.size());
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    if (x[k] == 0) {
      aux_mat = alpha.t() * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * ((exp(x[k] * beta) - 1.0) / beta)) * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}


//' Matrix-GEV density
//' 
//' Computes the density of a matrix-GEV distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' Does not allow for atoms in zero.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Transformation parameters.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector mgevden(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, Rcpp::NumericVector beta) {
  double mu{beta[0]};
  double sigma{beta[1]};
  double xi{beta[2]};
  
  Rcpp::NumericVector density(x.size());
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    if (xi == 0) {
      aux_mat = alpha.t() * matrix_exponential(S * exp(-(x[k] - mu) / sigma)) * exit_vect;
      density[k] = aux_mat(0,0) * exp(-(x[k] - mu) / sigma) / sigma;
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * pow(1.0 + (xi / sigma) * (x[k] - mu), -1.0 / xi)) * exit_vect;
      density[k] = aux_mat(0,0) * pow(1 + (xi / sigma) * (x[k] - mu), -(1 + xi) / xi) / sigma;
    }
  }
  return density;
}


//' Matrix-GEV cdf
//' 
//' Computes the cdf (tail) of a matrix-GEV distribution with parameters
//'  \code{alpha}, \code{S} and \code{beta} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Transformation parameters. 
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector mgevcdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, Rcpp::NumericVector beta, bool lower_tail = true) {
  double mu{beta[0]};
  double sigma{beta[1]};
  double xi{beta[2]};
  
  Rcpp::NumericVector cdf(x.size());
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    if (xi == 0) {
      aux_mat = alpha.t() * matrix_exponential(S * exp(-(x[k] - mu) / sigma)) * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * pow(1.0 + (xi / sigma) * (x[k] - mu), -1.0 / xi)) * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}


//' Discrete phase-type density
//' 
//' Computes the density of discrete phase-type distribution with parameters
//'  \code{alpha} and \code{S} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-transition matrix.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector dphdensity(Rcpp::NumericVector x, arma::vec alpha, arma::mat S) {
  Rcpp::NumericVector density(x.size());
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = e - (S * e);
  
  double max_val{max(x)};
  
  std::vector<arma::mat> vect = vector_of_powers(S, max_val);
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    aux_mat = alpha.t() * vect[x[k] - 1] * exit_vect;
    density[k] = aux_mat(0,0);
  }
  return density;
}


//' Discrete phase-type cdf
//' 
//' Computes the cdf (tail) of a discrete phase-type distribution with parameters 
//'  \code{alpha} and \code{S} at \code{x}.
//' 
//' @param x Non-negative value.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param lower_tail Cdf or tail.
//' @return The cdf (tail) at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector dphcdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, bool lower_tail = true) {
  Rcpp::NumericVector cdf(x.size());
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  
  double max_val{max(x)};
  
  std::vector<arma::mat> vect = vector_of_powers(S, max_val);
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < x.size(); ++k){
    aux_mat = alpha.t() * vect[x[k]] * e;
    cdf[k] = 1.0 - aux_mat(0,0);
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}

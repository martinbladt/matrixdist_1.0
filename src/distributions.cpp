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


//' Bivariate phase-type joint density of the feed forward type
//'
//' @param x Matrix of values.
//' @param alpha Vector of initial probabilities.
//' @param S11 Sub-intensity matrix.
//' @param S12 Matrix.
//' @param S22 Sub-intensity matrix.
//' @return Joint density at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector bivph_density(Rcpp::NumericMatrix x, arma::vec alpha, arma::mat S11, arma::mat S12, arma::mat S22) {
  long n{x.nrow()};
  
  Rcpp::NumericVector density(n);
  
  arma::mat e;
  e.ones(S22.n_cols, 1);
  arma::mat exit_vect = (S22 * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < n; ++k) {
    aux_mat = alpha.t() * matrix_exponential(S11 * x(k,0)) * S12 * matrix_exponential(S22 * x(k,1)) * exit_vect;
    density[k] = aux_mat(0,0);
  }
  return density;
}


//' Bivariate phase-type joint tail of the feed forward type
//'
//' @param x Matrix of values.
//' @param alpha Vector of initial probabilities.
//' @param S11 Sub-intensity matrix.
//' @param S12 Matrix.
//' @param S22 Sub-intensity matrix.
//' @return Joint tail at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector bivph_tail(Rcpp::NumericMatrix x, arma::vec alpha, arma::mat S11, arma::mat S12, arma::mat S22) {
  long n{x.nrow()};
  
  Rcpp::NumericVector tail(n);
  
  arma::mat e;
  e.ones(S22.n_cols, 1);
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < n; ++k) {
    aux_mat = alpha.t() * inv(S11 * (-1)) * matrix_exponential(S11 * x(k,0)) * S12 * matrix_exponential(S22 * x(k,1)) * e;
    tail[k] = aux_mat(0,0);
  }
  return tail;
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


//' Bivariate discrete phase-type joint density of the feed forward type
//'
//' @param x Matrix of values.
//' @param alpha Vector of initial probabilities.
//' @param S11 Sub-transition matrix.
//' @param S12 Matrix.
//' @param S22 Sub-transition matrix.
//' @return Joint density at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector bivdph_density(Rcpp::NumericMatrix x, arma::vec alpha, arma::mat S11, arma::mat S12, arma::mat S22) {
  long n{x.nrow()};
  
  Rcpp::NumericVector density(n);
  
  arma::mat e;
  e.ones(S22.n_cols, 1);
  arma::mat exit_vect = e - (S22 * e);
  
  double max_val1{max(x.column(0))};
  double max_val2{max(x.column(1))};
  
  std::vector<arma::mat> vect1 = vector_of_powers(S11, max_val1);
  std::vector<arma::mat> vect2 = vector_of_powers(S22, max_val2);
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < n; ++k) {
    aux_mat =  alpha.t() * vect1[x(k, 0) - 1] * S12 * vect2[x(k, 1) - 1] * exit_vect;
    density[k] = aux_mat(0,0);
  }
  return density;
}


//' Bivariate discrete phase-type joint tail of the feed forward type
//'
//' @param x Matrix of values.
//' @param alpha Vector of initial probabilities.
//' @param S11 Sub-transition matrix.
//' @param S12 Matrix.
//' @param S22 Sub-transition matrix.
//' @return Joint tail at \code{x}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector bivdph_tail(Rcpp::NumericMatrix x, arma::vec alpha, arma::mat S11, arma::mat S12, arma::mat S22) {
  long n{x.nrow()};
  
  Rcpp::NumericVector tail(n);
  
  arma::mat e;
  e.ones(S22.n_cols, 1);
  
  double max_val1{max(x.column(0))};
  double max_val2{max(x.column(1))};
  
  std::vector<arma::mat> vect1 = vector_of_powers(S11, max_val1);
  std::vector<arma::mat> vect2 = vector_of_powers(S22, max_val2);
  
  arma::mat aux_mat(1,1);
  
  for (int k{0}; k < n; ++k) {
    aux_mat =  alpha.t() * vect1[x(k, 0)] * S12 * vect2[x(k, 1)] * e;
    tail[k] = aux_mat(0,0);
  }
  return tail;
}


//' Multivariate discrete phase-type density
//' 
//' Computes the density of multivariate discrete phase-type distribution with 
//' parameters \code{alpha} and \code{S} at \code{x}.
//' 
//' @param x Matrix of positive integer values.
//' @param alpha Initial probabilities.
//' @param S_list List of marginal sub-transition matrices.
//' @return The density at \code{x}.
//' 
// [[Rcpp::export]]
arma::vec mdphdensity(Rcpp::NumericMatrix x, arma::vec alpha, Rcpp::List & S_list) {
  unsigned p{alpha.size()};
  long n{x.nrow()};
  long d{x.ncol()};
  
  arma::vec density(n);
  arma::mat aux_den(n, d);
  
  arma::mat e;
  e.ones(p, 1);
  
  std::vector<std::vector<arma::mat>> vect;
  std::vector<arma::mat> exit_vect; 
  
  for (int j{0}; j < d; ++j){
    double max_val{max(x.column(j))};
    arma::mat S = S_list[j];
    vect.push_back(vector_of_powers(S, max_val));
    exit_vect.push_back(e - (S * e));
  }
  
  
  arma::mat aux_mat(1,1);
  
  for (int i{0}; i < p; ++i){
    arma::mat in_vect(1, p);
    in_vect(0, i) = 1;
    for (int j{0}; j < d; ++j){
      for (int k{0}; k < n; ++k){
        aux_mat = in_vect * vect[j][x(k, j) - 1] * exit_vect[j];
        aux_den(k,j) = aux_mat(0,0);
      }
    }
    density = density + alpha[i] * arma::prod(aux_den, 1);
  }
  
  return density;
}

arma::mat kron_sum(arma::mat A, arma::mat B) {
  arma::mat I1;
  I1.eye(size(A));
  arma::mat I2;
  I2.eye(size(B));
  arma::mat K = kron(A, I2) + kron(I1, B);
  return K;
}


//' Bivariate phase-type joint density of the common shock type
//'
//' @param x Matrix of values.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-intensity matrix.
//' @param P Matrix.
//' @param Q1 Sub-intensity matrix.
//' @param Q2 Sub-intensity matrix.
//' @return Joint density at \code{x}.
//' @export
//' 
// [[Rcpp::export]]
Rcpp::NumericVector csph_density(Rcpp::NumericMatrix x, arma::vec alpha, arma::mat S, arma::mat P, arma::mat Q1, arma::mat Q2) {
  unsigned p1{S.n_rows};
  unsigned p2{Q1.n_rows};
  long n{x.nrow()};
  
  double m{0};
  
  Rcpp::NumericVector density(n);
  
  arma::mat e1;
  e1.ones(Q1.n_cols, 1);
  arma::mat exit_vect1 = (Q1 * (-1)) * e1;
  
  arma::mat e2;
  e2.ones(Q2.n_cols, 1);
  arma::mat exit_vect2 = (Q2 * (-1)) * e2;
  
  arma::mat exit_vect_prod = kron(exit_vect1, exit_vect2);
  arma::mat Q1pQ2 = kron_sum(Q1, Q2);
  
  arma::mat cmatrix(p2 * p2, p1);
  
  arma::mat aux_mat(1,1);
  
  std::vector<arma::mat> ei_vect;
  std::vector<arma::mat> eitei_vect; 
  
  for (int i{0}; i < p2; ++i) {
    arma::mat ei = arma::zeros(1,p2);
    ei[i] = 1;
    arma::mat eitei = kron(ei, ei);
    ei_vect.push_back(ei);
    eitei_vect.push_back(eitei);
  }
  
  
  for (int k{0}; k < n; ++k) {
    m = std::min(x(k,0), x(k,1));
    arma::mat B1 = matrix_exponential(Q1 * (x(k,0) - m));
    arma::mat B2 = matrix_exponential(Q2 * (x(k,1) - m));
    arma::mat B1tB2 = kron(B1, B2);
    arma::mat b_prod_alpha = B1tB2 * exit_vect_prod * alpha.t();
    
    arma::mat J = matrix_exponential(matrix_vanloan(Q1pQ2, S, b_prod_alpha) * m);
    for (int l{0}; l < p2 * p2; ++l) {
      for (int j{0}; j < p1; ++j) {
        cmatrix(l,j) = J(l,j + p2 * p2);
      }
    }
    
    for (int i{0}; i < p2; ++i) {
      aux_mat = eitei_vect[i] * cmatrix * P * ei_vect[i].t();
      density[k] += aux_mat(0,0);
    }
    
  }
  return density;
}



//' Bivariate phase-type joint density of the common shock type
//'
//' @param x Matrix of values.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-intensity matrix.
//' @param P Matrix.
//' @param Q1 Sub-intensity matrix.
//' @param Q2 Sub-intensity matrix.
//' @param h step size.
//' @return Joint density at \code{x}.
//' @export
//' 
// [[Rcpp::export]]
Rcpp::NumericVector csph_densityrk(Rcpp::NumericMatrix x, arma::vec alpha, arma::mat S, arma::mat P, arma::mat Q1, arma::mat Q2, double h) {
   unsigned p1{S.n_rows};
   unsigned p2{Q1.n_rows};
   long n{x.nrow()};
   
   double m{0};
   double eps{};
   
   Rcpp::NumericVector density(n);
   
   arma::mat e1;
   e1.ones(Q1.n_cols, 1);
   arma::mat exit_vect1 = (Q1 * (-1)) * e1;
   
   arma::mat e2;
   e2.ones(Q2.n_cols, 1);
   arma::mat exit_vect2 = (Q2 * (-1)) * e2;
   
   arma::mat k1(p2, p1);
   arma::mat k2(p2, p1);
   arma::mat k3(p2, p1);
   arma::mat k4(p2, p1);
   
   
   arma::mat aux_mat(1,1);
   
   for (int k{0}; k < n; ++k) {
     m = std::min(x(k,0), x(k,1));
     arma::mat B1_aux = matrix_exponential(Q1 * (x(k,0) - m)) * exit_vect1;
     arma::mat B2 = matrix_exponential(Q2 * (x(k,1) - m)) * exit_vect2 * alpha.t();
     
     int i{};
     i = m / h;
     
     double h2{};
     h2 = m / (i + 1);
     
     arma::mat mexp_inc1 = matrix_exponential(Q1 * (h2 / 2));
     arma::mat mexp_inc2 = matrix_exponential(Q2 * (h2 / 2));
     
     arma::mat mexp_cum1;
     arma::mat mexp_cum2;
     
     
     for (int i{0}; i < p2; ++i) {
       arma::mat ei = arma::zeros(1,p2);
       ei[i] = 1;
       arma::mat B1 = B1_aux * ei;
       arma::mat H = arma::zeros(p2, p1);
       mexp_cum1.eye(size(Q1));
       mexp_cum2.eye(size(Q2));
       for (eps = 0; eps <= m - h2 / 2; eps += h2) {
         k1 = H * S + mexp_cum1 * B1 * mexp_cum2 * B2;
         mexp_cum1 = mexp_cum1 * mexp_inc1;
         mexp_cum2 = mexp_cum2 * mexp_inc2;
         k2 = (H + k1 * (h2 / 2)) * S + mexp_cum1 * B1 * mexp_cum2 * B2;
         k3 = (H + k2 * (h2 / 2)) * S + mexp_cum1 * B1 * mexp_cum2 * B2;
         mexp_cum1 = mexp_cum1 * mexp_inc1;
         mexp_cum2 = mexp_cum2 * mexp_inc2;
         k4 = (H + k3 * h2) * S + mexp_cum1 * B1 * mexp_cum2 * B2;
         H = H  + (k1 + 2 * k2 + 2 * k3 + k4) * (h2 / 6);
       }
       
       aux_mat = ei * H * P * ei.t();
       density[k] += aux_mat(0,0);
     }
     
   }
   return density;
 }



//' Bivariate phase-type joint density of the common shock type
//'
//' @param x Matrix of values.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-intensity matrix.
//' @param P Matrix.
//' @param Q1 Sub-intensity matrix.
//' @param Q2 Sub-intensity matrix.
//' @return Joint density at \code{x}.
//' @export
//' 
// [[Rcpp::export]]
Rcpp::NumericVector csph_density_test(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, arma::mat P, arma::mat Q1, arma::mat Q2) {
  unsigned p1{S.n_rows};
  unsigned p2{Q1.n_rows};
   
  double m{0};
   
  Rcpp::NumericVector density(1);
   
  arma::mat e1;
  e1.ones(Q1.n_cols, 1);
  arma::mat exit_vect1 = (Q1 * (-1)) * e1;
   
  arma::mat e2;
  e2.ones(Q2.n_cols, 1);
  arma::mat exit_vect2 = (Q2 * (-1)) * e2;
   
  arma::mat exit_vect_prod = kron(exit_vect1, exit_vect2);
  arma::mat Q1pQ2 = kron_sum(Q1, Q2);
   
  arma::mat cmatrix(p2 * p2, p1);
   
  arma::mat aux_mat(1,1);
   
  m = std::min(x[0], x[1]);
  arma::mat B1 = matrix_exponential(Q1 * (x[0] - m));
  arma::mat B2 = matrix_exponential(Q2 * (x[1] - m));
  arma::mat B1tB2 = kron(B1, B2);
  arma::mat b_prod_alpha = B1tB2 * exit_vect_prod * alpha.t();
  
  arma::mat J = matrix_exponential(matrix_vanloan(Q1pQ2, S, b_prod_alpha) * m);
  for (int l{0}; l < p2 * p2; ++l) {
    for (int j{0}; j < p1; ++j) {
      cmatrix(l,j) = J(l,j + p2 * p2);
    }
  }
  
  for (int i{0}; i < p2; ++i) {
     arma::mat ei = arma::zeros(1,p2);
     ei[i] = 1;
     arma::mat eitei = kron(ei, ei);
     aux_mat = eitei * cmatrix * P * ei.t();
     density[0] += aux_mat(0,0);
  }
     
  return density;
 }


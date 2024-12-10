#ifndef EM_LL_UNI_PI   // if x.h hasn't been included yet...
#define EM_LL_UNI_PI   //   #define this so the compiler knows it has been included

#include <RcppArmadillo.h>

double logLikelihood_UNIs_PI(double h, arma::vec & alpha,
                             arma::mat & S,
                             SEXP beta1,
                             SEXP beta2,
                             const Rcpp::NumericVector & obs,
                             const Rcpp::NumericVector & weight,
                             SEXP rcens, const Rcpp::NumericVector & rcweight,
                             const Rcpp::NumericVector & scale1,
                             const Rcpp::NumericVector & scale2,
                             const std::string & gfun_name);
  
#endif
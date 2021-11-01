# include <RcppArmadillo.h>
#include "exp_arm.h"
// [[ Rcpp :: depends ( RcppArmadillo )]]

//' Product of two matrices
//' 
//' @param A1 Matrix.
//' @param A2 Matrix.
//' @return Computes C = A1 * A2.
//' 
// [[Rcpp::export]]
Rcpp::NumericMatrix matrix_product(Rcpp::NumericMatrix A1, Rcpp::NumericMatrix A2) {
  arma::mat AA1 = Rcpp::as<arma::mat>(A1);
  arma::mat AA2 = Rcpp::as<arma::mat>(A2);
  return(Rcpp::wrap(AA1 * AA2));
}

//' Inverse of a matrix
//' 
//' @param A A matrix.
//' 
// [[Rcpp::export]]
Rcpp::NumericMatrix matrix_inverse(Rcpp::NumericMatrix A) {
  arma::mat AA = Rcpp::as<arma::mat>(A);
  return(Rcpp::wrap(inv(AA)));
}

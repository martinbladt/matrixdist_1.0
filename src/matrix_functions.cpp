#include <Rcpp.h>
using namespace Rcpp;


//' Clone a vector 
//' 
//' @param v A vector.
//' 
// [[Rcpp::export]]
NumericVector clone_vector(NumericVector v) {
  NumericVector new_v = clone(v);
  return new_v;
}

//' Clone a matrix 
//' 
//' @param m A matrix.
//' 
// [[Rcpp::export]]
NumericMatrix clone_matrix(NumericMatrix m) {
  NumericMatrix new_m = clone(m);
  return new_m;
}





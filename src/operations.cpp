#include <Rcpp.h>
using namespace Rcpp;
#include "matrix_functions.h"

//' Computes the initial distribution and sub-intensity of the sum of PH
//' 
//' @param pi1 initial distribution
//' @param T1 sub-intensity
//' @param pi2 initial distribution
//' @param T2 sub-intensity
//' 
// [[Rcpp::export]]
List sumPH(NumericVector pi1, NumericMatrix T1, NumericVector pi2, NumericMatrix T2) {
  
  long p1{T1.nrow()};
  long p2{T2.nrow()};
  
  NumericVector pi_sum(p1 + p2);
  NumericMatrix T_sum(p1 + p2, p1 + p2);
  
  NumericMatrix pi2_m(1, p2, pi2.begin());
  
  NumericVector m_e(p1, 1);
  NumericMatrix e(p1, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T1 * (-1), e);
  
  NumericMatrix auxmat = matrix_product(t, pi2_m);
  
  
  for (int i{0}; i < p1 + p2; ++i) {
    if (i < p1) {
      pi_sum[i] = pi1[i];
    }
    
    for (int j{0}; j < p1 + p2; ++j) {
      if (i < p1) {
        if (j < p1) {
          T_sum(i,j) = T1(i,j);
        }
        else {
          T_sum(i,j) = auxmat(i,j - p1);
        }
      }
      else if (i >= p1 && j>= p1) {
        T_sum(i,j) = T2(i - p1,j - p1);
      }
    }
  }
  
  List L = List::create(Named("pi") = pi_sum, _["T"] = T_sum);
  
  return L;
}

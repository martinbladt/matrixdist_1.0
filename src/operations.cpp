#include <Rcpp.h>
using namespace Rcpp;
#include "matrix_functions.h"

//' Computes the initial distribution and sub-intensity of the sum of PH
//' 
//' @param alpha1 initial distribution
//' @param S1 sub-intensity
//' @param alpha2 initial distribution
//' @param S2 sub-intensity
//' 
// [[Rcpp::export]]
List sumPH(NumericVector alpha1, NumericMatrix S1, NumericVector alpha2, NumericMatrix S2) {
  
  long p1{S1.nrow()};
  long p2{S2.nrow()};
  
  NumericVector alpha_sum(p1 + p2);
  NumericMatrix S_sum(p1 + p2, p1 + p2);
  
  NumericMatrix alpha2_m(1, p2, alpha2.begin());
  
  NumericVector m_e(p1, 1);
  NumericMatrix e(p1, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S1 * (-1), e);
  
  NumericMatrix auxmat = matrix_product(t, alpha2_m);
  
  
  for (int i{0}; i < p1 + p2; ++i) {
    if (i < p1) {
      alpha_sum[i] = alpha1[i];
    }
    
    for (int j{0}; j < p1 + p2; ++j) {
      if (i < p1) {
        if (j < p1) {
          S_sum(i,j) = S1(i,j);
        }
        else {
          S_sum(i,j) = auxmat(i,j - p1);
        }
      }
      else if (i >= p1 && j>= p1) {
        S_sum(i,j) = S2(i - p1,j - p1);
      }
    }
  }
  
  List L = List::create(Named("alpha") = alpha_sum, _["S"] = S_sum);
  
  return L;
}

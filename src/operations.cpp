#include <Rcpp.h>
using namespace Rcpp;
#include "matrix_functions.h"


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


// [[Rcpp::export]]
NumericMatrix Kroneckerproduct(NumericMatrix a, NumericMatrix b) 
{ 
  int arows{a.nrow()};
  int acolumns{a.nrow()};
  int brows{b.nrow()};
  int bcolumns{b.ncol()};
  NumericMatrix c(arows * brows, acolumns * bcolumns);
  for (size_t i = 0; i < arows; ++i)
    for (size_t j = 0; j < acolumns; ++j)
      for (size_t k = 0; k < brows; ++k)
        for (size_t l = 0; l < bcolumns; ++l)
          c(i*brows + k, j*bcolumns + l) = a(i, j) * b(k, l);
  return c;
} 


// [[Rcpp::export]]
NumericMatrix Kroneckersum(NumericMatrix a, NumericMatrix b) 
{ 
  return matrix_sum(Kroneckerproduct(a, NumericMatrix::diag(b.nrow(), 1.0)), Kroneckerproduct(NumericMatrix::diag(a.nrow(), 1.0), b));
}

#include "auxilliary.h"
# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]



// [[Rcpp::export]]
double LInf_normArma(arma::mat A) {
  double value{0.0};
  
  for (int i{0}; i < A.n_rows; ++i) {
    double row_sum{0.0};
    for (int j{0}; j < A.n_cols; ++j) {
      row_sum += abs(A(i,j));
    }
    value = std::max(value, row_sum);
  }
  return value;
}



//' Creates the matrix  (A1, B1 ; 0, A2)
//' @param A1 matrix
//' @param A2 matrix
//' @param B1 matrix
//' @return Computes (A1, B1 ; 0, A2)
//' 
// [[Rcpp::export]]
arma::mat matrix_VanLoanArma(arma::mat A1, arma::mat A2, arma::mat B1) {
  unsigned p1{A1.n_rows};
  unsigned p2{A2.n_rows};
  unsigned p{p1 + p2};
  
  arma::mat auxiliarMatrix(p, p);
  
  for (int i{0}; i < p; ++i) {
    for (int j{0}; j < p; ++j) {
      if ( i < p1 && j < p1) {
        auxiliarMatrix(i,j) = A1(i,j);
      }
      else if (i >= p1 && j < p1) {
        auxiliarMatrix(i,j) = 0;
      }
      else if (i < p1 && j >= p1) {
        auxiliarMatrix(i,j) = B1(i,j - p1);
      }
      else {
        auxiliarMatrix(i,j) = A2(i - p1,j - p1);
      }
    }
  }
  return auxiliarMatrix;
}


// [[Rcpp::export]]
double matrixMaxDiagonal_arma(const arma::mat & A) {
  double maximum{A(0,0)};
  for (int i{0}; i < A.n_rows; ++i) {
    if (A(i,i) > maximum) {
      maximum = A(i,i);
    }
  }
  return maximum;
}

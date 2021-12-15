#include "m_exp.h"
# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]


//' L inf norm of a matrix
//' 
//' Computes the L inf norm of a matrix \code{A}, which is defined as:
//' L_inf A =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
//' 
//' @param A A matrix.
//' @return The L inf norm.
//' 
// [[Rcpp::export]]
double inf_norm(arma::mat A) {
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
//' 
//' @param A1 Matrix.
//' @param A2 Matrix.
//' @param B1 Matrix.
//' @return Computes (A1, B1 ; 0, A2).
//' 
// [[Rcpp::export]]
arma::mat matrix_VanLoan(arma::mat A1, arma::mat A2, arma::mat B1) {
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


//' Maximum diagonal element of a matrix
//' 
//' @param A Matrix.
//' @return The maximum value in the diagonal.
//' 
// [[Rcpp::export]]
double max_diagonal(const arma::mat & A) {
  double maximum{A(0,0)};
  for (int i{0}; i < A.n_rows; ++i) {
    if (A(i,i) > maximum) {
      maximum = A(i,i);
    }
  }
  return maximum;
}


//' Matrix exponential algorithm
//' 
//' MATLAB's built-in algorithm - Pade approximation.
//' 
//' @param A A matrix.
//' 
// [[Rcpp::export]]
arma::mat matrix_exponential(arma::mat A) {
  const int q{6};
  
  arma::mat matrixAuxiliar(A.n_rows,A.n_cols);
  arma::mat ExpM(A.n_rows,A.n_cols);
  
  double aNorm{inf_norm(A)};
  
  int ee{static_cast<int>(log2(aNorm)) + 1};
  
  int s{std::max(0, ee + 1)};
  
  double t{1.0 / pow(2.0, s)};
  
  arma::mat a2 = A * t;
  arma::mat x = a2;
  
  double c{0.5};
  
  ExpM.eye( size(A) );
  
  ExpM = ExpM + (a2 * c);
  
  arma::mat d; d.eye( size(A) );
  
  d = (d + a2 * (-c));
  
  int p{1};
  
  for (int k{2}; k <= q; ++k) {
    c = c * static_cast<double>(q - k + 1) / static_cast<double>(k * (2 * q - k + 1));
    
    x = (a2 * x);
    
    ExpM = (x * c) + ExpM;
    
    if (p) {
      d = (x * c) + d;
    }
    else {
      d = (x * (-c)) + d;
    }
    p = !p;
  }
  //  E -> inverse(D) * E
  ExpM = inv(d) * ExpM;
  //  E -> E^(2*S)
  for (int k = 1; k <= s; ++k) {
    ExpM = ExpM * ExpM;
  }
  return(ExpM);
}


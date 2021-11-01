#include "matrix_functions.h"
#include "exp_arm.h"
#include <Rcpp.h>
using namespace Rcpp;

// The main objective is to implement a matrix exponential method and other functions for matrix that we need afterwards


//' Add matrices
//' 
//' Computes C =  A + B 
//' @param A A matrix
//' @param B A matrix
//' 
// [[Rcpp::export]]
NumericMatrix matrix_sum(const NumericMatrix & A, const NumericMatrix & B) {
  long rows = A.nrow();
  long cols = A.ncol();
  
  NumericMatrix C(rows, cols);
  for (int i{0}; i < rows; ++i) {
    for (int j{0}; j < cols; ++j) {
      C(i,j) = A(i,j) + B(i,j);
    }
  }
  return (C);
}

//' L-oo norm of a matrix
//' 
//' Computes the L-oo norm of a matrix \code{A}, which is defined as:
//' L-oo A =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
//' @param A A matrix
// [[Rcpp::export]]
double LInf_norm(const NumericMatrix & A) {
  double value{0.0};
  
  for (int i{0}; i < A.nrow(); ++i) {
    double row_sum{0.0};
    for (int j{0}; j < A.ncol(); ++j) {
      row_sum += std::abs(A(i,j));
    }
    value = std::max(value, row_sum);
  }
  return value;
}


//' Maximum entry in a matrix
//' 
//' Find the maximum entry
//' @param A a matrix
//' 
// [[Rcpp::export]]
double matrixMax(const NumericMatrix & A) {
  double maximum{A(0,0)};
  for (int i{0}; i < A.nrow(); ++i) {
    for (int j{0}; j < A.ncol(); ++j) {
      if (A(i,j) > maximum) {
        maximum = A(i,j);
      }
    }
  }
  return maximum;
}

//' Maximum entry in the diagonal of a matrix
//' 
//' @param A a matrix
// [[Rcpp::export]]
double matrixMaxDiagonal(const NumericMatrix & A) {
  double maximum{A(0,0)};
  for (int i{0}; i < A.nrow(); ++i) {
    if (A(i,i) > maximum) {
      maximum = A(i,i);
    }
  }
  return maximum;
}

//' Computes A^n
//' 
//' @param n integer
//' @param A a matrix
//' 
// [[Rcpp::export]]
NumericMatrix matrix_power(int n, const NumericMatrix & A) {
  if (n == 1) {
    return A;
  }
  else if (n == 2) {
    return matrix_product(A, A);
  }
  NumericMatrix previousMatrix(matrix_product(A, A));
  NumericMatrix newMatrix(matrix_product(A, previousMatrix));
  for (int i{4}; i <= n; ++i) {
    previousMatrix = clone(newMatrix);
    newMatrix = matrix_product(A, previousMatrix);
  }
  return newMatrix;
}

//' Clone a vector 
//' 
//' @param v a vector
//' 
// [[Rcpp::export]]
NumericVector clone_vector(NumericVector v) {
  NumericVector new_v = clone(v);
  return new_v;
}

//' Clone a matrix 
//' 
//' @param m a matrix
//' 
// [[Rcpp::export]]
NumericMatrix clone_matrix(NumericMatrix m) {
  NumericMatrix new_m = clone(m);
  return new_m;
}


//' Creates the matrix  (A1, B1 ; 0, A2)
//' 
//' @param A1 a matrix
//' @param A2 a matrix
//' @param B1 a matrix
//' 
// [[Rcpp::export]]
NumericMatrix matrix_VanLoan(const NumericMatrix & A1, const NumericMatrix & A2, const NumericMatrix & B1) {
  long p1{A1.nrow()};
  long p2{A2.nrow()};
  long p{p1 + p2};
  
  NumericMatrix auxiliarMatrix(p, p);
  
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

//' Creates a matrix with the given vector in the diagonal
//' 
//' @param vec a vector
//' 
// [[Rcpp::export]]
NumericMatrix diagonal_vector(const NumericVector & vec) {
  NumericMatrix diagonalMatrix(vec.size(),vec.size());
  for (int i{0}; i < vec.size(); ++i) {
    diagonalMatrix(i,i) = vec[i];
  }
  return diagonalMatrix;
}


#include "matrix_functions.h"
#include <Rcpp.h>
using namespace Rcpp;

// The main objective is to implement a matrix exponential method and other functions for matrix that we need afterwards

//' Product of two matrices
//' @param a matrix
//' @param b matrix
//' @return Computes c = a * b
//' 
// [[Rcpp::export]]
NumericMatrix matrix_product(NumericMatrix a, NumericMatrix b) {
  
  const size_t rows = a.nrow(), cols = b.ncol(), n = a.ncol();
  
  NumericMatrix c(rows, cols);
  
  for (size_t i = 0; i < rows; ++i) {
    auto ci = c(i,_);
    for (size_t j = 0; j < cols; ++j) {
      ci[j] = 0;
    }
  }
  
  for (size_t i = 0; i < rows; ++i) {
    const auto ai = a(i,_);
    auto ci = c(i,_);
    
    for (size_t k = 0; k < n; ++k) {
      const auto bk = b(k,_);
      const auto aik = ai[k];
      
      for (size_t j = 0; j < cols; ++j) {
        ci[j] += aik * bk[j];
      }
    }
  }
  return (c);
}

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
      row_sum += abs(A(i,j));
    }
    value = std::max(value, row_sum);
  }
  return value;
}

//' Solves a system with multiple right hand sides
//' 
//' AX=B which can be decomposed as LUX=B and finds X.
//' When B is the identity matrix the solution is the inverse of A
//' @param A1 a matrix
//' @param B a matrix
//' 
// [[Rcpp::export]]
NumericMatrix solve_linear_system(NumericMatrix A1, const NumericMatrix & B) {
  long i{};
  long ipiv{};
  long j{};
  long jcol{};
  double piv{};
  double t{};
  
  NumericMatrix A = clone(A1);
  NumericMatrix X = clone(B);
  
  for (jcol = 1; jcol <= A.nrow(); ++jcol) {
    //  Find the maximum element in column I.
    piv = abs(A(jcol - 1,jcol - 1));
    ipiv = jcol;
    for (i = jcol + 1; i <= A.nrow(); ++i) {
      if (piv < abs(A(i - 1,jcol - 1))) {
        piv = abs(A(i - 1,jcol - 1));
        ipiv = i;
      }
    }
    
    if (piv == 0.0) {
      Rcerr << "\n";
      Rcerr << "R8MAT_FSS_NEW - Fatal error!\n";
      Rcerr << "  Zero pivot on step " << jcol << "\n";
      Rcpp::stop("Error!\n"); 
    }
    
    //  Switch rows JCOL and IPIV, and X.
    if ( jcol != ipiv ) {
      for ( j = 1; j <= A.nrow(); ++j ) {
        t              = A(jcol-1, j-1);
        A(jcol-1, j-1) = A(ipiv-1, j-1);
        A(ipiv-1, j-1) = t;
      }
      for ( j = 0; j < A.nrow(); ++j ) {
        t            = X(jcol-1, j);
        X(jcol-1, j) = X(ipiv-1, j);
        X(ipiv-1, j) = t;
      }
    }
    
    //  Scale the pivot row.
    t = A(jcol - 1,jcol - 1);
    A(jcol - 1,jcol - 1) = 1.0;
    for (j = jcol + 1; j <= A.nrow(); ++j) {
      A(jcol - 1,j - 1) /= t;
    }
    for (j = 0; j < B.ncol(); ++j) {
      X(jcol - 1,j) /= t;
    }
    //  Use the pivot row to eliminate lower entries in that column.
    for (i = jcol + 1; i <= A.nrow(); ++i) {
      if (A(i - 1,jcol - 1) != 0.0) {
        t = - A(i - 1,jcol - 1);
        A(i - 1,jcol - 1) = 0.0;
        for (j = jcol + 1; j <= A.nrow(); ++j) {
          A(i - 1,j - 1) +=  t * A(jcol - 1,j - 1);
        }
        for (j = 0; j < B.ncol(); ++j) {
          X(i - 1,j) +=  t * X(jcol - 1,j);
        }
      }
    }
  }
  
  //  Back solve.
  for (jcol = A.nrow(); 2 <= jcol; --jcol) {
    for (i = 1; i < jcol; ++i) {
      for (j = 0; j < B.ncol(); ++j) {
        X(i - 1,j) -= A(i - 1,jcol - 1) * X(jcol - 1,j);
      }
    }
  }
  return X;
}

//' Inverse of a matrix
//' 
//' Computes the inverse
//' @param A a matrix
//' 
// [[Rcpp::export]]
NumericMatrix matrix_inverse(NumericMatrix A) {
  return solve_linear_system(A, NumericMatrix::diag(A.nrow(), 1.0));
}


//' Matrix exponential algorithm
//' 
//' MATLAB's built-in algorithm - Pade approximation
//' @param A a matrix
//' 
// [[Rcpp::export]]
NumericMatrix matrix_exponential(const NumericMatrix & A) {
  
  const int q{6};
  
  NumericMatrix matrixAuxiliar(A.nrow(),A.ncol());
  NumericMatrix ExpM(A.nrow(),A.ncol());
  
  double aNorm{LInf_norm(A)};
  
  int ee{static_cast<int>(log2(aNorm)) + 1};
  
  int s{std::max(0, ee + 1)};
  
  double t{1.0 / pow(2.0, s)};
  
  
  NumericMatrix a2 = clone(A * t);
  NumericMatrix x = clone(a2);
  
  
  double c{0.5};
  
  ExpM = NumericMatrix::diag(A.nrow(), 1.0);
  
  ExpM = matrix_sum(ExpM, a2 * c);
  
  NumericMatrix d = NumericMatrix::diag(A.nrow(), 1.0);
  
  d = matrix_sum(d, a2 * (-c));
  
  int p{1};
  
  for (int k{2}; k <= q; ++k) {
    c = c * static_cast<double>(q - k + 1) / static_cast<double>(k * (2 * q - k + 1));
    
    x = matrix_product(a2, x);
    
    ExpM = matrix_sum(x * c, ExpM);
    
    if (p) {
      d = matrix_sum(x * c, d);
    }
    else {
      d = matrix_sum(x * (-c), d);
    }
    p = !p;
  }
  //  E -> inverse(D) * E
  ExpM = solve_linear_system(d, ExpM);
  //  E -> E^(2*S)
  for (int k = 1; k <= s; ++k) {
    ExpM = matrix_product(ExpM, ExpM);
  }
  return ExpM;
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


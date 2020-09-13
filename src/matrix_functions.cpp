#include <Rcpp.h>
using namespace Rcpp;

// The main objective is to implement a matrix exponential method and other functions for matrix that we need afterwards


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

//' Computes C =  A + B 
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

//' Returns the  L-oo norm of a matrix
//' 
//' The matrix L-oo norm is defined as:
//' L-oo A =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
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
//' AX=B which can be decompose as LUX=B and finds X
//' When B is the identity matrix the solution is the inverse of A
// [[Rcpp::export]]
NumericMatrix solve_linear_system(NumericMatrix A, const NumericMatrix & B ) {
  long i{};
  long ipiv{};
  long j{};
  long jcol{};
  double piv{};
  double t{};
  
  NumericMatrix X = B;
  
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
      std::cerr << "\n";
      std::cerr << "R8MAT_FSS_NEW - Fatal error!\n";
      std::cerr << "  Zero pivot on step " << jcol << "\n";
      exit (1);
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
    for (j = jcol+1; j <= A.nrow(); ++j) {
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


// [[Rcpp::export]]
NumericMatrix matrix_inverse(const NumericMatrix & A) {
  return solve_linear_system(A, NumericMatrix::diag(A.nrow(), 1.0));
}


//' MATLAB's built-in matrix exponential algorithm
// [[Rcpp::export]]
NumericMatrix matrix_exponential(const NumericMatrix & A ) {
  
  const int q{6};
  
  NumericMatrix matrixAuxiliar(A.nrow(),A.ncol());
  NumericMatrix ExpM(A.nrow(),A.ncol());
  
  double aNorm{LInf_norm(A)};
  
  int ee{static_cast<int>(log2(aNorm)) + 1};
  
  int s{std::max(0, ee + 1)};
  
  double t{1.0 / pow(2.0, s)};
  
  
  NumericMatrix a2 = A * t;
  NumericMatrix x = a2;
  
  
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


// [[Rcpp::export]]
double matrixMax(const NumericMatrix & A) {
  // Returns the maximum value in a matrix
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

// [[Rcpp::export]]
double matrixMaxDiagonal(const NumericMatrix & A) {
  // Returns the maximum value in the diagonal of a matrix
  double maximum{A(0,0)};
  for (int i{0}; i < A.nrow(); ++i) {
    if (A(i,i) > maximum) {
      maximum = A(i,i);
    }
  }
  return maximum;
}

//' Computes A^n
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
    previousMatrix = newMatrix;
    newMatrix = matrix_product(A, previousMatrix);
  }
  return newMatrix;
}


// [[Rcpp::export]]
NumericVector phdensity(NumericVector x, NumericVector pi, NumericMatrix T) {
  
  NumericVector density(x.size());
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  NumericVector e(pi.size(), 1);
  NumericMatrix m_e(pi.size(), 1, e.begin());
  NumericMatrix m_t = matrix_product(T * (-1), m_e);
  
  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      density[k] = (1.0 - matrix_product(m_pi, m_e)(0,0));
    }
    else {
      density[k] = (matrix_product(m_pi, matrix_product(matrix_exponential(T * x[k]), m_t))(0,0));
    }
  }
  return density;
}


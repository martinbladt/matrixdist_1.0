#include "auxilliary.h"
# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]


// [[Rcpp::export]]
Rcpp::NumericMatrix productArma(Rcpp::NumericMatrix A1, Rcpp::NumericMatrix A2) {
  arma::mat AA1 = Rcpp::as<arma::mat>(A1);
  arma::mat AA2 = Rcpp::as<arma::mat>(A2);
  return(Rcpp::wrap(AA1 * AA2));
}


// [[Rcpp::export]]
Rcpp::NumericMatrix sumArma(arma::mat A1, arma::mat A2) {
  return(Rcpp::wrap(A1 + A2));
}


// [[Rcpp::export]]
Rcpp::NumericMatrix invArma(arma::mat A) {
  return(Rcpp::wrap(inv(A)));
}


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

//' Matrix exponential algorithm
//' 
//' MATLAB's built-in algorithm - Pade approximation
//' @param Ainput a matrix
//' 
// [[Rcpp::export]]
arma::mat mexponentialArma(arma::mat Ainput) {
  
  arma::mat A = (Ainput); 
  
  const int q{6};
  
  arma::mat matrixAuxiliar(A.n_rows,A.n_cols);
  arma::mat ExpM(A.n_rows,A.n_cols);
  
  double aNorm{LInf_normArma(A)};
  
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

//' Creates the matrix  (A1, B1 ; 0, A2)
//' @param A1 matrix
//' @param A2 matrix
//' @param B1 matrix
//' @return Computes (A1, B1 ; 0, A2)
//' 
// [[Rcpp::export]]
arma::mat matrix_VanLoanArma(arma::mat A1, arma::mat A2, arma::mat B1) {
  long p1{A1.n_rows};
  long p2{A2.n_rows};
  long p{p1 + p2};
  
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

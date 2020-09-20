#include <Rcpp.h>
using namespace Rcpp;
#include "matrix_functions.h"

//' Random structure of a phase-type
//' 
//' Generates random parameters \code{pi} and \code{T} of a phase-type distribution of dimension \code{p} with chosen structure
//' @param p Dimension of the phase-type
//' @param structure Type of structure: "General", "Hyperexponential", "GErlang", "Coxian" or "GCoxian"
//' @param scale_factor A factor that multiplies the sub-intensity matrix
//' @return Random parameters \code{pi} and \code{T} of a phase-type
//' @examples
//' random_structure(3) 
//' random_structure(5, "Hyperexponential") 
// [[Rcpp::export]]
List random_structure(int p, String structure = "General", double scale_factor = 1) {
  // Structure of pi and T
  NumericVector pi_legal(p);
  NumericMatrix T_legal(p, p);
  
  if (structure == "General") {
    for (int i = 0; i < p; ++i) {
      pi_legal[i] = 1;
      for (int j = 0; j < p; ++j) {
        T_legal(i, j) = 1;
      }
    }
  }
  else if (structure == "Hyperexponential") {
    for (int i = 0; i < p; i++) {  
      pi_legal[i] = 1;
      T_legal(i, i) = 1;
    }
  }
  else if (structure == "GErlang") {
    pi_legal[0] = 1;
    for (int i = 0; i < p - 1; ++i) {
      T_legal(i, i + 1) = 1;
    } 
    T_legal(p - 1, p - 1) = 1;
  }
  else if (structure == "Coxian") {
    pi_legal[0] = 1;
    for (int i = 0; i < p - 1; ++i) { 
      T_legal(i, i) = 1;
      T_legal(i, i + 1) = 1;
    }
    T_legal(p - 1, p - 1) = 1; 
  }
  else if (structure == "GCoxian") {
    for (int i = 0; i < p - 1; ++i) { 
      pi_legal[i] = 1;
      T_legal(i, i) = 1;
      T_legal(i, i + 1) = 1;
    }
    T_legal(p - 1, p - 1) = 1;
    pi_legal[p - 1] = 1;
  }
  else{
    Rcerr << "Non-existent structure\n";
  }
  
  
  // Random initialization
  NumericVector pi_rnd(p);
  NumericMatrix T_rnd(p, p);
  
  double sum{0.0};
  
  for (int i{0}; i < p; ++i) {
    if (pi_legal[i] == 1) {
      pi_rnd[i] = runif(1)[0];
      sum += pi_rnd[i];
    }
  }
  
  for (int i{0}; i < p; ++i) {
    pi_rnd[i] = pi_rnd[i] / sum;
  }
  
  for (int i{0}; i < p; ++i) {
    for (int j{0}; j < p; ++j) {
      if ((i != j) && (T_legal(i,j) == 1)) {
        T_rnd(i,j) = runif(1)[0];
        T_rnd(i,i) -= T_rnd(i,j);
      }
    }
  }
  
  // exit rate
  for (int i{0}; i < p; ++i) {
    if (T_legal(i,i) == 1) {
      T_rnd(i,i) -= runif(1)[0];
    }
  }
  T_rnd = T_rnd * scale_factor;
  
  List L = List::create(Named("pi") = pi_rnd , _["T"] = T_rnd);
  
  return L;
}

//' Random reward matrix
//' 
//' The rows of the matrix sum to 1
// [[Rcpp::export]]
NumericMatrix random_reward(int p, int dim) {
  NumericMatrix R(p,dim);
  
  for (int i{0}; i < p; ++i) {
    double sum = 0;
    for (int j{0}; j < dim; ++j) {
        R(i,j) = runif(1)[0];
        sum += R(i,j);
    }
    for (int j{0}; j < dim; ++j) {
      R(i,j) = R(i,j) / sum;
    }
  }
  return R;
}


//' Changes the values of T and R to make the MPH* having rewards that sum to 1
//' 
//' Better clone T and R?
// [[Rcpp::export]]
void norm_mph(NumericMatrix T, NumericMatrix R){
  NumericVector sumrow(T.nrow());
  for (int i{0}; i < R.nrow(); ++i) {
    for (int j{0}; j < R.ncol(); ++j) {
      sumrow[i] = sumrow[i] + R(i,j);
    }
  }
  sumrow = 1.0 / sumrow;
  
  NumericMatrix diagMat(diagonal_vector(sumrow));
  
  T = matrix_product(diagMat, T);
  R = matrix_product(diagMat, R);
  
}



//' Merges the matrices T11, T12 and T22 of a bivariate matrix
//' 
// [[Rcpp::export]]
NumericMatrix merge_matrices(NumericMatrix T11, NumericMatrix T12, NumericMatrix T22) {
  long p1{T11.nrow()};
  long p2{T22.nrow()};
  NumericMatrix T(p1 +p2, p1 + p2);

  for (int i{0}; i < p1; ++i) {
    for (int j{0}; j < p1; ++j) {
      T(i,j) = T11(i,j);
    }
    for (int j{0}; j < p2; ++j) {
      T(i,j + p1) = T12(i,j);
    }
  }
  for (int i{0}; i < p2; ++i) {
    for (int j{0}; j < p2; ++j) {
      T(i + p1,j + p1) = T22(i,j);
    }
  }
  
  return T;
}


//' Random parameters of a bivariate phase-type
// [[Rcpp::export]]
List random_phase_BivPH(int p1, int p2, double scale_factor = 1) {
  
  NumericVector alpha(p1);
  NumericMatrix T11(p1, p1);
  NumericMatrix T12(p1, p2);
  NumericMatrix T22(p2, p2);
  
  NumericMatrix R(p1 + p2, 2);
  
  double sum{0.0};
  
  for (int i{0}; i < p1; ++i) {
    alpha[i] = runif(1)[0];
    sum += alpha[i];
  }
  alpha = alpha / sum;
  
  for (int i{0}; i < p1; ++i) {
    for (int j{0}; j < p1; ++j) {
      if ((i != j)) {
        T11(i,j) = runif(1)[0];
        T11(i,i) -= T11(i,j);
      }
    }
    for (int j{0}; j < p2; ++j) {
      T12(i,j) = runif(1)[0];
      T11(i,i) -= T12(i,j);
    }
    R(i,0) = 1;
  }
  
  for (int i{0}; i < p2; ++i) {
    for (int j{0}; j < p2; ++j) {
      if ((i != j)) {
        T22(i,j) = runif(1)[0];
        T22(i,i) -= T22(i,j);
      }
    }
    R(i + p1,1) = 1;
  }
  
  for (int i{0}; i < p2; ++i) {
    T22(i,i) -= runif(1)[0];
  }
  
  T11 = T11 * scale_factor;
  T12 = T12 * scale_factor;
  T22 = T22 * scale_factor;
  
  NumericMatrix T(merge_matrices(T11, T12, T22));
  
  NumericVector pi(p1 + p2);
  for (int i{0}; i < p1; ++i) {
    pi[i] = alpha[i];
  }
  
  List L = List::create(Named("alpha") = alpha, _["T11"] = T11, _["T12"] = T12, _["T22"] = T22, _["pi"] = pi, _["T"] = T, _["R"] = R);
  
  return L;
}


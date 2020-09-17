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


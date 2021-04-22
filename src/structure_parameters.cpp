#include <Rcpp.h>
using namespace Rcpp;
#include "matrix_functions.h"
#include "exp_arm.h"

//' Random structure of a phase-type
//' 
//' Generates random parameters \code{alpha} and \code{S} of a phase-type distribution of dimension \code{p} with chosen structure
//' @param p Dimension of the phase-type
//' @param structure Type of structure: "general", "hyperexponential", "gerlang", "coxian" or "gcoxian"
//' @param scale_factor A factor that multiplies the sub-intensity matrix
//' @return Random parameters \code{alpha} and \code{S} of a phase-type
//' 
// [[Rcpp::export]]
List random_structure(int p, String structure = "general", double scale_factor = 1) {
  // Structure of alpha and S
  NumericVector alpha_legal(p);
  NumericMatrix S_legal(p, p);
  
  if (structure == "general") {
    for (int i = 0; i < p; ++i) {
      alpha_legal[i] = 1;
      for (int j = 0; j < p; ++j) {
        S_legal(i, j) = 1;
      }
    }
  }
  else if (structure == "hyperexponential") {
    for (int i = 0; i < p; i++) {  
      alpha_legal[i] = 1;
      S_legal(i, i) = 1;
    }
  }
  else if (structure == "gerlang") {
    alpha_legal[0] = 1;
    for (int i = 0; i < p - 1; ++i) {
      S_legal(i, i + 1) = 1;
    } 
    S_legal(p - 1, p - 1) = 1;
  }
  else if (structure == "coxian") {
    alpha_legal[0] = 1;
    for (int i = 0; i < p - 1; ++i) { 
      S_legal(i, i) = 1;
      S_legal(i, i + 1) = 1;
    }
    S_legal(p - 1, p - 1) = 1; 
  }
  else if (structure == "gcoxian") {
    for (int i = 0; i < p - 1; ++i) { 
      alpha_legal[i] = 1;
      S_legal(i, i) = 1;
      S_legal(i, i + 1) = 1;
    }
    S_legal(p - 1, p - 1) = 1;
    alpha_legal[p - 1] = 1;
  }
  else{
    Rcerr << "non-existent structure\n";
  }
  
  
  // Random initialization
  NumericVector alpha_rnd(p);
  NumericMatrix S_rnd(p, p);
  
  double sum{0.0};
  
  for (int i{0}; i < p; ++i) {
    if (alpha_legal[i] == 1) {
      alpha_rnd[i] = runif(1)[0];
      sum += alpha_rnd[i];
    }
  }
  
  for (int i{0}; i < p; ++i) {
    alpha_rnd[i] = alpha_rnd[i] / sum;
  }
  
  for (int i{0}; i < p; ++i) {
    for (int j{0}; j < p; ++j) {
      if ((i != j) && (S_legal(i,j) == 1)) {
        S_rnd(i,j) = runif(1)[0];
        S_rnd(i,i) -= S_rnd(i,j);
      }
    }
  }
  
  // exit rate
  for (int i{0}; i < p; ++i) {
    if (S_legal(i,i) == 1) {
      S_rnd(i,i) -= runif(1)[0];
    }
  }
  S_rnd = S_rnd * scale_factor;
  
  List L = List::create(Named("alpha") = alpha_rnd , _["S"] = S_rnd);
  
  return L;
}

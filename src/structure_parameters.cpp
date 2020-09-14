#include <Rcpp.h>
using namespace Rcpp;

// Posible structures: General, Hyperexponential, GErlang, Coxian, GCoxian

// [[Rcpp::export]]
NumericVector random_structure(int p, String structure = "General") {
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
  return pi_legal;
}




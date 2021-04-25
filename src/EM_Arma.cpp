#include "auxilliary.h"
# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]

////////////////////////////////////////////
// Uniformization
////////////////////////////////////////////

//' Computes elements S^n / n! until the value size
//' @param theVector a vector
//' @param S sub-untensity matrix
//' @param a a number
//' @param sizevect size of vector
// [[Rcpp::export]]
void vectorOfMatrices_arma(std::vector<arma::mat> & theVector, const arma::mat & S, double a, int sizevect) {
  long p{S.n_rows};
  
  arma::mat I;
  I.eye(size(S));
  
  arma::mat P = I + S *(1/ a);
  
  theVector.push_back(I);
  
  for (int k{1}; k <= sizevect; ++k) {
    theVector.push_back( (P * (1.0 / k) ) * theVector[k - 1]);
  }
}


//' Computes e^(Sx) base on the values on powerVector
//' @param x a number
//' @param n an integer
//' @param powerVector a vector
//' @param a a number
// [[Rcpp::export]]
arma::mat matrixExpSum_arma(double x, int n, const std::vector<arma::mat> & powerVector, double a) {
  arma::mat resultmatrix = powerVector[0];
  
  for (int i{1}; i <= n; ++i) {
    resultmatrix = resultmatrix + powerVector[i] * exp(i * log(a * x));
  }
  resultmatrix = resultmatrix * exp(-a * x);
  
  return resultmatrix;
}

// Computes (A)^(2^n)
// [[Rcpp::export]]
void pow2Matrix_arma(int n , arma::mat & A) {
  arma::mat auxMat(A.n_rows,A.n_cols);
  
  for (int i{1}; i <= n; ++i) {
    auxMat = A * A;
    A = auxMat;
  }
}

//  Find n such that P(N > n) = epsilon with N Poisson distributed
int findN_arma(double epsilon, double lambda) {
  int n{0};
  double cumProb{0.0};
  
  do {
    cumProb += R::dpois(n, lambda, false);
    ++n;
  } while (cumProb < 1.0 - epsilon);
  
  return (n - 1);
}


//' EM using Uniformization for matrix exponential
//' 
//' @param epsilon positive parameter
//' @param alpha initial probalities
//' @param S sub-intensity
//' @param obs the observations
//' @param weight the weights for the observations
//' @param rcens censored observations
//' @param rcweight the weights for the censored observations
//' 
// [[Rcpp::export]]
void EMstep_UNI(double epsilon, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  long p{S.n_rows};
  
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat Bmean = arma::zeros(p,1);
  arma::mat Zmean = arma::zeros(p,1);
  arma::mat Nmean = arma::zeros(p,p + 1);
  
  arma::mat avector(1,p);
  arma::mat bvector(p,1);
  arma::mat cmatrix(p,p);
  arma::mat aux_exp(p,p);
  
  arma::mat aux_mat(1,1);
  
  arma::mat J(2 * p,2 * p);
  arma::mat tProductPi(p,p);
  tProductPi = t * alpha.t();
  
  J = matrix_VanLoanArma(S, S, tProductPi);
  
  double a = matrixMaxDiagonal_arma(J * (-1));
  
  int N{findN_arma(epsilon, 1)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma(theVector, J, a, N);
  
  
  double SumOfWeights{0.0};
  double density{0.0};
  
  //E-step
  //  Unccensored data
  for (int k{0}; k < obs.size(); ++k) {
    SumOfWeights += weight[k];
    
    double x{obs[k]};
    
    if (x * a <= 1.0) {
      J = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = log(a * x) / log(2);
      ++n;
      
      J = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, J);
    }
    
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }
    
    avector = alpha.t() * aux_exp;
    bvector = aux_exp * t;
    aux_mat = alpha.t() * bvector;
    density = aux_mat(0,0);
    
    //E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha[i] * bvector(i,0) * weight[k] / density;
      Nmean(i,p) += avector(0,i) * t(i,0) * weight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * weight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += S(i,j) * cmatrix(j,i) * weight[k] / density;
      }
    }
  }
  //  Right-Censored Data
  double SumOfCensored{0.0};
  if (rcens.size() > 0) {
    tProductPi = e * alpha.t();
    J = matrix_VanLoanArma(S, S, tProductPi);
    theVector.clear();
    vectorOfMatrices_arma(theVector, J, a, N);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    SumOfCensored += rcweight[k];
    
    double x{rcweight[k]};
    
    if (x * a <= 1.0) {
      J = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = log(a * x) / log(2);
      ++n;
      
      J = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, J);
    }
    
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }
    
    bvector = aux_exp * e;
    aux_mat = alpha.t() * bvector;
    density = aux_mat(0,0);
    
    //E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha[i] * bvector(i,0) * rcweight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * rcweight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += S(i,j) * cmatrix(j,i) * rcweight[k] / density;
      }
    }
  }
  
  // M step
  for (int i{0}; i < p; ++i) {
    alpha[i] = Bmean(i,0) / (SumOfWeights + SumOfCensored);
    if (alpha[i] < 0) {
      alpha[i] = 0;
    }
    t(i,0) = Nmean(i,p) / Zmean(i,0);
    if (t(i,0) < 0) {
      t(i,0) = 0;
    }
    S(i,i) = -t(i,0);
    for (int j{0}; j < p; ++j) {
      if (i != j) {
        S(i,j) = Nmean(i,j) / Zmean(i,0);
        if (S(i,j) < 0) {
          S(i,j) = 0;
        }
        S(i,i) -= S(i,j);
      }
    }
  }
}

////////////////////////////////////////////
// Modified EM with Matlab algorithm
////////////////////////////////////////////

//' Computes elements S^n / n! until the value size
//' @param theVector a vector
//' @param S sub-untensity matrix
//' @param sizevect size of vector
// [[Rcpp::export]]
void vectorOfMatrices_arma2(std::vector<arma::mat> & theVector, const arma::mat & S, int sizevect) {
  
  arma::mat I;
  I.eye(size(S));
  
  theVector.push_back(I);
  
  for (int k{1}; k <= sizevect; ++k) {
    theVector.push_back( S * theVector[k - 1]);
  }
}



//' EM using Matlab algorithm for matrix exponential in combination with Armadillo
//' 
//' @param epsilon nuisance parameter
//' @param alpha initial probalities
//' @param S sub-intensity
//' @param obs the observations
//' @param weight the weights for the observations
//' @param rcens censored observations
//' @param rcweight the weights for the censored observations
//' 
// [[Rcpp::export]]
void EMstep_PADE(double epsilon, arma::vec & alpha,  arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  long p{S.n_rows};
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat Bmean = arma::zeros(p,1);
  arma::mat Zmean = arma::zeros(p,1);
  arma::mat Nmean = arma::zeros(p,p + 1);
  
  arma::mat avector(1,p);
  arma::mat bvector(p,1);
  arma::mat cmatrix(p,p);
  arma::mat aux_exp(p,p);
  
  arma::mat aux_mat(1,1);
  
  arma::mat J(2 * p,2 * p);
  arma::mat tProductPi(p,p);
  tProductPi = t * alpha.t();
  
  
  J = matrix_VanLoanArma(S, S, tProductPi);
  
  
  double JNorm{LInf_normArma(J)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma2(theVector, J, 6);
  
  
  arma::mat X(2 * p,2 * p);
  arma::mat D(2 * p,2 * p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  
  
  double SumOfWeights{0.0};
  double density{0.0};
  
  //E-step
  //  Unccensored data
  for (int k{0}; k < obs.size(); ++k) {
    SumOfWeights += weight[k];
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * obs[k])) + 1};
    s = std::max(0, ee + 1);
    xmod = obs[k] / pow(2, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    J = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      J = J + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    J = inv(D) * J;
    for (int l = 1; l <= s; ++l) {
      J = J * J;
    }
    
    // Separate matrix
    
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }
    
    avector = alpha.t() * aux_exp;
    bvector = aux_exp * t;
    aux_mat = alpha.t() * bvector;
    density = aux_mat(0,0);
    
    //E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha[i] * bvector(i,0) * weight[k] / density;
      Nmean(i,p) += avector(0,i) * t(i,0) * weight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * weight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += S(i,j) * cmatrix(j,i) * weight[k] / density;
      }
    }
  }
  //  Right-Censored Data
  double SumOfCensored{0.0};
  if (rcens.size() > 0) {
    tProductPi = e * alpha.t();
  }
  for (int k{0}; k < rcens.size(); ++k) {
    SumOfCensored += rcweight[k];
    
    J = mexponentialArma(matrix_VanLoanArma(S, S, tProductPi) * rcens[k]);
    
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }
    
    bvector = aux_exp * e;
    aux_mat = alpha.t() * bvector;
    density = aux_mat(0,0);
    
    //E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha[i] * bvector(i,0) * rcweight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * rcweight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += S(i,j) * cmatrix(j,i) * rcweight[k] / density;
      }
    }
  }
  
  // M step
  for (int i{0}; i < p; ++i) {
    alpha[i] = Bmean(i,0) / (SumOfWeights + SumOfCensored);
    if (alpha[i] < 0) {
      alpha[i] = 0;
    }
    t(i,0) = Nmean(i,p) / Zmean(i,0);
    if (t(i,0) < 0) {
      t(i,0) = 0;
    }
    S(i,i) = -t(i,0);
    for (int j{0}; j < p; ++j) {
      if (i != j) {
        S(i,j) = Nmean(i,j) / Zmean(i,0);
        if (S(i,j) < 0) {
          S(i,j) = 0;
        }
        S(i,i) -= S(i,j);
      }
    }
  }
}

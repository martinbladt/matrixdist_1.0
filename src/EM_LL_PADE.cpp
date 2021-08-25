#include "auxilliary.h"
# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]


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
//' @param h nuisance parameter
//' @param alpha initial probalities
//' @param S sub-intensity
//' @param obs the observations
//' @param weight the weights for the observations
//' @param rcens censored observations
//' @param rcweight the weights for the censored observations
//' 
// [[Rcpp::export]]
void EMstep_PADE(double h, arma::vec & alpha,  arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  unsigned p{S.n_rows};
  
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
    xmod = obs[k] / pow(2.0, s);
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
    J = matrix_VanLoanArma(S, S, tProductPi);
    JNorm = LInf_normArma(J);
    theVector.clear();
    vectorOfMatrices_arma2(theVector, J, 6);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    SumOfCensored += rcweight[k];
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * rcens[k])) + 1};
    s = std::max(0, ee + 1);
    xmod = rcens[k] / pow(2.0, s);
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
// Loglikelihoods
////////////////////////////////////////////


//' Loglikelihood of PH using Pade
//' 
//' Loglikelihood for a sample 
//' @param h nuisance parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodPH_PADE(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  unsigned p{S.n_rows};
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat E(size(S));
  
  double JNorm{LInf_normArma(S)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma2(theVector, S, 6);
  
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * obs[k])) + 1};
    s = std::max(0, ee + 1);
    xmod = obs[k] / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * t;
    density = aux_mat(0,0);
    logLh += weight[k] * std::log(density);
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * rcens[k])) + 1};
    s = std::max(0, ee + 1);
    xmod = rcens[k] / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-Weibull using Pade
//' 
//' Loglikelihood for a sample 
//' @param h nuisance parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta in-homogeneity parameter
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMweibull_PADE(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat E(size(S));
  
  double JNorm{LInf_normArma(S)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma2(theVector, S, 6);
  
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * pow(obs[k], beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = pow(obs[k], beta) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta) + (beta - 1) * std::log(obs[k]));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * pow(rcens[k], beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = pow(rcens[k], beta) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}




//' Loglikelihood of matrix-Pareto using Pade
//' 
//' Loglikelihood for a sample 
//' @param h nuisance parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta in-homogeneity parameter
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMpareto_PADE(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat E(size(S));
  
  double JNorm{LInf_normArma(S)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma2(theVector, S, 6);
  
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * std::log(obs[k] / beta + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = std::log(obs[k] / beta + 1) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) - std::log(obs[k] + beta));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * std::log(rcens[k] / beta + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = std::log(rcens[k] / beta + 1) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}




//' Loglikelihood of matrix-lognormal using Pade
//' 
//' Loglikelihood for a sample 
//' @param h nuisance parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta in-homogeneity parameter
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMlognormal_PADE(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat E(size(S));
  
  double JNorm{LInf_normArma(S)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma2(theVector, S, 6);
  
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * pow(std::log(obs[k] + 1), beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = pow(std::log(obs[k] + 1), beta) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta) + (beta -1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * pow(std::log(rcens[k] + 1), beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = pow(std::log(rcens[k] + 1), beta) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-loglogistic using Pade
//' 
//' Loglikelihood for a sample 
//' @param h nuisance parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta in-homogeneity parameter
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMloglogistic_PADE(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta[0] < 0 || beta[1] < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat E(size(S));
  
  double JNorm{LInf_normArma(S)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma2(theVector, S, 6);
  
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * std::log(pow(obs[k] / beta[0], beta[1]) + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = std::log(pow(obs[k] / beta[0], beta[1]) + 1) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta[1]) - std::log(beta[0]) + (beta[1] - 1) * (std::log(obs[k]) - std::log(beta[0])) - std::log(pow(obs[k] / beta[0], beta[1]) + 1));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * std::log(pow(rcens[k] / beta[0], beta[1]) + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = std::log(pow(rcens[k] / beta[0], beta[1]) + 1) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-Gompertz using Pade
//' 
//' Loglikelihood for a sample 
//' @param h nuisance parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta in-homogeneity parameter
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMgompertz_PADE(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat E(size(S));
  
  double JNorm{LInf_normArma(S)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma2(theVector, S, 6);
  
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * (exp(obs[k] * beta) - 1) / beta)) + 1};
    s = std::max(0, ee + 1);
    xmod = (exp(obs[k] * beta) - 1) / beta / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + obs[k] * beta);
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * (exp(rcens[k] * beta) - 1) / beta)) + 1};
    s = std::max(0, ee + 1);
    xmod = (exp(rcens[k] * beta) - 1) / beta / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}

//' Loglikelihood of matrix-GEV using Pade
//' 
//' Loglikelihood for a sample 
//' @param h nuisance parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta in-homogeneity parameter
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMgev_PADE(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta[1] < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat E(size(S));
  
  double JNorm{LInf_normArma(S)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma2(theVector, S, 6);
  
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  if (beta[2] == 0) {
    // Non censored data
    for (int k{0}; k < obs.size(); ++k) {
      
      // Matrix exponential
      int pind{1};
      int ee{static_cast<int>(log2(JNorm  * exp(-(obs[k] - beta[0]) / beta[1]))) + 1};
      s = std::max(0, ee + 1);
      xmod = exp(-(obs[k] - beta[0]) / beta[1])/ pow(2.0, s);
      c = 0.5;
      X = theVector[1] * (c * xmod);
      
      E = theVector[0] + X;
      D = theVector[0] - X;
      
      for (int l{2}; l <= q; ++l) {
        c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
        X = theVector[l] * (c * pow(xmod,l));
        E = E + X;
        if (pind) {
          D =  D + X;
        }
        else {
          D = D - X;
        }
        pind = !pind;
      }
      E = inv(D) * E;
      for (int l = 1; l <= s; ++l) {
        E = E * E;
      }
      
      aux_mat = alpha.t() * E * t;
      density = aux_mat(0,0);
      logLh += weight[k] * (std::log(density) - std::log(beta[1]) - (obs[k] - beta[0]) / beta[1]);
    }
    //Right censored data
    for (int k{0}; k < rcens.size(); ++k) {
      
      // Matrix exponential
      int pind{1};
      int ee{static_cast<int>(log2(JNorm  * exp(-(rcens[k] - beta[0]) / beta[1]))) + 1};
      s = std::max(0, ee + 1);
      xmod = exp(-(rcens[k] - beta[0]) / beta[1]) / pow(2.0, s);
      c = 0.5;
      X = theVector[1] * (c * xmod);
      
      E = theVector[0] + X;
      D = theVector[0] - X;
      
      for (int l{2}; l <= q; ++l) {
        c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
        X = theVector[l] * (c * pow(xmod,l));
        E = E + X;
        if (pind) {
          D =  D + X;
        }
        else {
          D = D - X;
        }
        pind = !pind;
      }
      E = inv(D) * E;
      for (int l = 1; l <= s; ++l) {
        E = E * E;
      }
      
      aux_mat = alpha.t() * E * e;
      density = aux_mat(0,0);
      logLh += rcweight[k] * std::log(density);
    }
  }
  else{
    // Non censored data
    for (int k{0}; k < obs.size(); ++k) {
      
      // Matrix exponential
      int pind{1};
      int ee{static_cast<int>(log2(JNorm  * pow(1 + (beta[2] / beta[1]) * (obs[k] - beta[0]) , - 1 / beta[2]))) + 1};
      s = std::max(0, ee + 1);
      xmod = pow(1 + (beta[2] / beta[1]) * (obs[k] - beta[0]) , - 1 / beta[2]) / pow(2.0, s);
      c = 0.5;
      X = theVector[1] * (c * xmod);
      
      E = theVector[0] + X;
      D = theVector[0] - X;
      
      for (int l{2}; l <= q; ++l) {
        c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
        X = theVector[l] * (c * pow(xmod,l));
        E = E + X;
        if (pind) {
          D =  D + X;
        }
        else {
          D = D - X;
        }
        pind = !pind;
      }
      E = inv(D) * E;
      for (int l = 1; l <= s; ++l) {
        E = E * E;
      }
      
      aux_mat = alpha.t() * E * t;
      density = aux_mat(0,0);
      logLh += weight[k] * (std::log(density) - std::log(beta[1]) - (1 + 1 / beta[2]) * std::log(1 + (beta[2] / beta[1]) * (obs[k] - beta[0])));
    }
    //Right censored data
    for (int k{0}; k < rcens.size(); ++k) {
      
      // Matrix exponential
      int pind{1};
      int ee{static_cast<int>(log2(JNorm  * pow(1 + (beta[2] / beta[1]) * (rcens[k] - beta[0]) , - 1 / beta[2]))) + 1};
      s = std::max(0, ee + 1);
      xmod = pow(1 + (beta[2] / beta[1]) * (rcens[k] - beta[0]) , - 1 / beta[2]) / pow(2.0, s);
      c = 0.5;
      X = theVector[1] * (c * xmod);
      
      E = theVector[0] + X;
      D = theVector[0] - X;
      
      for (int l{2}; l <= q; ++l) {
        c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
        X = theVector[l] * (c * pow(xmod,l));
        E = E + X;
        if (pind) {
          D =  D + X;
        }
        else {
          D = D - X;
        }
        pind = !pind;
      }
      E = inv(D) * E;
      for (int l = 1; l <= s; ++l) {
        E = E * E;
      }
      
      aux_mat = alpha.t() * E * e;
      density = aux_mat(0,0);
      logLh += rcweight[k] * std::log(density);
    }
  }
  return logLh;
}


////////////////////////////////////////////
// Scaled versions of loglikelihoods (for regression):
////////////////////////////////////////////


//' Loglikelihood of PH using Pade
//' 
//' Loglikelihood for a sample 
//' @param h nuisance parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' @param scale1 scale for observations
//' @param scale2 scale for censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodPH_PADEs(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  unsigned p{S.n_rows};
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat E(size(S));
  
  double JNorm{LInf_normArma(S)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma2(theVector, S, 6);
  
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * scale1[k] * obs[k])) + 1};
    s = std::max(0, ee + 1);
    xmod = scale1[k] * obs[k] / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * scale2[k] * rcens[k])) + 1};
    s = std::max(0, ee + 1);
    xmod = scale2[k] * rcens[k] / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-Weibull using Pade
//' 
//' Loglikelihood for a sample 
//' @param h nuisance parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta in-homogeneity parameter
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' @param scale1 scale for observations
//' @param scale2 scale for censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMweibull_PADEs(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat E(size(S));
  
  double JNorm{LInf_normArma(S)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma2(theVector, S, 6);
  
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * scale1[k] * pow(obs[k], beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale1[k] * pow(obs[k], beta) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta) + (beta - 1) * std::log(obs[k]));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * scale2[k] * pow(rcens[k], beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale2[k] * pow(rcens[k], beta) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}




//' Loglikelihood of matrix-Pareto using Pade
//' 
//' Loglikelihood for a sample 
//' @param h nuisance parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta in-homogeneity parameter
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' @param scale1 scale for observations
//' @param scale2 scale for censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMpareto_PADEs(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat E(size(S));
  
  double JNorm{LInf_normArma(S)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma2(theVector, S, 6);
  
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * scale1[k] * std::log(obs[k] / beta + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale1[k] * std::log(obs[k] / beta + 1) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) - std::log(obs[k] + beta));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * scale2[k] * std::log(rcens[k] / beta + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale2[k] * std::log(rcens[k] / beta + 1) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-lognormal using Pade
//' 
//' Loglikelihood for a sample 
//' @param h nuisance parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta in-homogeneity parameter
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' @param scale1 scale for observations
//' @param scale2 scale for censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMlognormal_PADEs(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat E(size(S));
  
  double JNorm{LInf_normArma(S)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma2(theVector, S, 6);
  
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * scale1[k] * pow(std::log(obs[k] + 1), beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale1[k] * pow(std::log(obs[k] + 1), beta) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta) + (beta -1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * scale2[k] * pow(std::log(rcens[k] + 1), beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale2[k] * pow(std::log(rcens[k] + 1), beta) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-loglogistic using Pade
//' 
//' Loglikelihood for a sample 
//' @param h nuisance parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta in-homogeneity parameter
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' @param scale1 scale for observations
//' @param scale2 scale for censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMloglogistic_PADEs(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta[0] < 0 || beta[1] < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat E(size(S));
  
  double JNorm{LInf_normArma(S)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma2(theVector, S, 6);
  
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * scale1[k] * std::log(pow(obs[k] / beta[0], beta[1]) + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale1[k] * std::log(pow(obs[k] / beta[0], beta[1]) + 1) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta[1]) - std::log(beta[0]) + (beta[1] - 1) * (std::log(obs[k]) - std::log(beta[0])) - std::log(pow(obs[k] / beta[0], beta[1]) + 1));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * scale2[k] * std::log(pow(rcens[k] / beta[0], beta[1]) + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale2[k] * std::log(pow(rcens[k] / beta[0], beta[1]) + 1) / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-Gompertz using Pade
//' 
//' Loglikelihood for a sample 
//' @param h nuisance parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta in-homogeneity parameter
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' @param scale1 scale for observations
//' @param scale2 scale for censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMgompertz_PADEs(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat E(size(S));
  
  double JNorm{LInf_normArma(S)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma2(theVector, S, 6);
  
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * scale1[k] * (exp(obs[k] * beta) - 1) / beta)) + 1};
    s = std::max(0, ee + 1);
    xmod = scale1[k] * (exp(obs[k] * beta) - 1) / beta / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + obs[k] * beta);
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * scale2[k] * (exp(rcens[k] * beta) - 1) / beta)) + 1};
    s = std::max(0, ee + 1);
    xmod = scale2[k] * (exp(rcens[k] * beta) - 1) / beta / pow(2.0, s);
    c = 0.5;
    X = theVector[1] * (c * xmod);
    
    E = theVector[0] + X;
    D = theVector[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = theVector[l] * (c * pow(xmod,l));
      E = E + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    E = inv(D) * E;
    for (int l = 1; l <= s; ++l) {
      E = E * E;
    }
    
    aux_mat = alpha.t() * E * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


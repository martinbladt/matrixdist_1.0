#include "auxilliary.h"
# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]


////////////////////////////////////////////
// EM Uniformization 
////////////////////////////////////////////

//' Computes elements S^n / n! until the value size
//' @param theVector a vector
//' @param S sub-untensity matrix
//' @param a a number
//' @param sizevect size of vector
// [[Rcpp::export]]
void vectorOfMatrices_arma(std::vector<arma::mat> & theVector, const arma::mat & S, double a, int sizevect) {
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
    resultmatrix = resultmatrix + powerVector[i] * exp(i * std::log(a * x));
  }
  resultmatrix = resultmatrix * exp(-a * x);
  
  return resultmatrix;
}

// Computes (A)^(2^n)
// [[Rcpp::export]]
void pow2Matrix_arma(int n , arma::mat & A) {
  arma::mat auxMat(size(A));
  
  for (int i{1}; i <= n; ++i) {
    auxMat = A * A;
    A = auxMat;
  }
}

//  Find n such that P(N > n) = h with N Poisson distributed
int findN_arma(double h, double lambda) {
  int n{0};
  double cumProb{0.0};
  
  do {
    cumProb += R::dpois(n, lambda, false);
    ++n;
  } while (cumProb < 1.0 - h);
  
  return (n - 1);
}


//' EM using Uniformization for matrix exponential
//' 
//' @param h positive parameter
//' @param alpha initial probalities
//' @param S sub-intensity
//' @param obs the observations
//' @param weight the weights for the observations
//' @param rcens censored observations
//' @param rcweight the weights for the censored observations
//' 
// [[Rcpp::export]]
void EMstep_UNI(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
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
  
  double a = matrixMaxDiagonal_arma(J * (-1));
  
  int N{findN_arma(h, 1)};
  
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
      n = std::log(a * x) / std::log(2.0);
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
    
    double x{rcens[k]};
    
    if (x * a <= 1.0) {
      J = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
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
// Log-likelihoods
////////////////////////////////////////////



//' Loglikelihood using Uniformization
//' 
//' Loglikelihood for a sample 
//' @param h positive parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodPH_UNI(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat s = (S * (-1)) * e;
  
  arma::mat mExp(size(S));
  
  double a = matrixMaxDiagonal_arma(S * (-1));
  
  int N{findN_arma(h, 1)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma(theVector, S, a, N);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{obs[k]};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * s;
    density = aux_mat(0,0);
    logLh += weight[k] * std::log(density);
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{rcens[k]};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of matrix-Weibull using Uniformization
//' 
//' Loglikelihood for a sample 
//' @param h positive parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMweibull_UNI(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  
  if(beta < 0) return NA_REAL;
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat s = (S * (-1)) * e;
  
  arma::mat mExp(size(S));
  
  double a = matrixMaxDiagonal_arma(S * (-1));
  
  int N{findN_arma(h, 1)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma(theVector, S, a, N);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{pow(obs[k], beta)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * s;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta) + (beta - 1) * std::log(obs[k]));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{pow(rcens[k], beta)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-Pareto using Uniformization
//' 
//' Loglikelihood for a sample 
//' @param h positive parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMpareto_UNI(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  
  if(beta < 0) return NA_REAL;
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat s = (S * (-1)) * e;
  
  arma::mat mExp(size(S));
  
  double a = matrixMaxDiagonal_arma(S * (-1));
  
  int N{findN_arma(h, 1)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma(theVector, S, a, N);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{std::log(obs[k] / beta + 1)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * s;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) - std::log(obs[k] + beta));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{std::log(rcens[k] / beta + 1)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-lognormal using Uniformization
//' 
//' Loglikelihood for a sample 
//' @param h positive parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMlognormal_UNI(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  
  if(beta < 0) return NA_REAL;
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat s = (S * (-1)) * e;
  
  arma::mat mExp(size(S));
  
  double a = matrixMaxDiagonal_arma(S * (-1));
  
  int N{findN_arma(h, 1)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma(theVector, S, a, N);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{pow(std::log(obs[k] + 1), beta)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * s;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta) + (beta -1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{pow(std::log(rcens[k] + 1), beta)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-loglogistic using Uniformization
//' 
//' Loglikelihood for a sample 
//' @param h positive parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMloglogistic_UNI(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  
  if(beta[0] < 0 || beta[1] < 0) return NA_REAL;
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat s = (S * (-1)) * e;
  
  arma::mat mExp(size(S));
  
  double a = matrixMaxDiagonal_arma(S * (-1));
  
  int N{findN_arma(h, 1)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma(theVector, S, a, N);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{std::log(pow(obs[k] / beta[0], beta[1]) + 1)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * s;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta[1]) - std::log(beta[0]) + (beta[1] - 1) * (std::log(obs[k]) - std::log(beta[0])) - std::log(pow(obs[k] / beta[0], beta[1]) + 1));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{std::log(pow(rcens[k] / beta[0], beta[1]) + 1)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-Gompertz using Uniformization
//' 
//' Loglikelihood for a sample 
//' @param h positive parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMgompertz_UNI(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  
  if(beta < 0) return NA_REAL;
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat s = (S * (-1)) * e;
  
  arma::mat mExp(size(S));
  
  double a = matrixMaxDiagonal_arma(S * (-1));
  
  int N{findN_arma(h, 1)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma(theVector, S, a, N);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{(exp(obs[k] * beta) - 1) / beta};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * s;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + obs[k] * beta);
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{(exp(rcens[k] * beta) - 1) / beta};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-GEV using Uniformization
//' 
//' Loglikelihood for a sample 
//' @param h positive parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMgev_UNI(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  
  if(beta[1] < 0) return NA_REAL;
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat s = (S * (-1)) * e;
  
  arma::mat mExp(size(S));
  
  double a = matrixMaxDiagonal_arma(S * (-1));
  
  int N{findN_arma(h, 1)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma(theVector, S, a, N);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  if (beta[2] == 0) {
    // Non censored data
    for (int k{0}; k < obs.size(); ++k) {
      double x{exp(-(obs[k] - beta[0]) / beta[1])};
      
      if (x * a <= 1.0) {
        mExp = matrixExpSum_arma(x, N, theVector, a);
      }
      else {
        int n{};
        n = std::log(a * x) / std::log(2.0);
        ++n;
        
        mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
        
        pow2Matrix_arma(n, mExp);
      }
      aux_mat = alpha.t() * mExp * s;
      density = aux_mat(0,0);
      logLh += weight[k] * (std::log(density) - std::log(beta[1]) - (obs[k] - beta[0]) / beta[1]);
    }
    //Right censored data
    for (int k{0}; k < rcens.size(); ++k) {
      double x{exp(-(rcens[k] - beta[0]) / beta[1])};
      
      if (x * a <= 1.0) {
        mExp = matrixExpSum_arma(x, N, theVector, a);
      }
      else {
        int n{};
        n = std::log(a * x) / std::log(2.0);
        ++n;
        
        mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
        
        pow2Matrix_arma(n, mExp);
      }
      aux_mat = alpha.t() * mExp * e;
      density = aux_mat(0,0);
      logLh += rcweight[k] * std::log(density);
    }
  }
  else {
    // Non censored data
    for (int k{0}; k < obs.size(); ++k) {
      double x{pow(1 + (beta[2] / beta[1]) * (obs[k] - beta[0]) , - 1 / beta[2])};
      
      if (x * a <= 1.0) {
        mExp = matrixExpSum_arma(x, N, theVector, a);
      }
      else {
        int n{};
        n = std::log(a * x) / std::log(2.0);
        ++n;
        
        mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
        
        pow2Matrix_arma(n, mExp);
      }
      aux_mat = alpha.t() * mExp * s;
      density = aux_mat(0,0);
      logLh += weight[k] * (std::log(density) - std::log(beta[1]) - (1 + 1 / beta[2]) * std::log(1 + (beta[2] / beta[1]) * (obs[k] - beta[0])));
    }
    //Right censored data
    for (int k{0}; k < rcens.size(); ++k) {
      double x{pow(1 + (beta[2] / beta[1]) * (rcens[k] - beta[0]) , - 1 / beta[2])};
      
      if (x * a <= 1.0) {
        mExp = matrixExpSum_arma(x, N, theVector, a);
      }
      else {
        int n{};
        n = std::log(a * x) / std::log(2.0);
        ++n;
        
        mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
        
        pow2Matrix_arma(n, mExp);
      }
      aux_mat = alpha.t() * mExp * e;
      density = aux_mat(0,0);
      logLh += rcweight[k] * std::log(density);
    }
  }
  return logLh;
}


////////////////////////////////////////////
// Scaled versions of loglikelihoods (for regression):
////////////////////////////////////////////


//' Loglikelihood of PH using Uniformization
//' 
//' Loglikelihood for a sample 
//' @param h positive parameter
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
double logLikelihoodPH_UNIs(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat s = (S * (-1)) * e;
  
  arma::mat mExp(size(S));
  
  double a = matrixMaxDiagonal_arma(S * (-1));
  
  int N{findN_arma(h, 1)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma(theVector, S, a, N);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{scale1[k] * obs[k]};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * s;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{scale2[k] * rcens[k]};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-Weibull using Uniformization
//' 
//' Loglikelihood for a sample 
//' @param h positive parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' @param scale1 scale for observations
//' @param scale2 scale for censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMweibull_UNIs(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  
  if(beta < 0) return NA_REAL;
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat s = (S * (-1)) * e;
  
  arma::mat mExp(size(S));
  
  double a = matrixMaxDiagonal_arma(S * (-1));
  
  int N{findN_arma(h, 1)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma(theVector, S, a, N);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{scale1[k] * pow(obs[k], beta)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * s;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta) + (beta -1) * std::log(obs[k]));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{scale2[k] * pow(rcens[k], beta)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-Pareto using Uniformization
//' 
//' Loglikelihood for a sample 
//' @param h positive parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' @param scale1 scale for observations
//' @param scale2 scale for censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMpareto_UNIs(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  
  if(beta < 0) return NA_REAL;
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat s = (S * (-1)) * e;
  
  arma::mat mExp(size(S));
  
  double a = matrixMaxDiagonal_arma(S * (-1));
  
  int N{findN_arma(h, 1)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma(theVector, S, a, N);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{scale1[k] * std::log(obs[k] / beta + 1)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * s;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) - std::log(obs[k] + beta));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{scale2[k] * std::log(rcens[k] / beta + 1)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-lognormal using Uniformization
//' 
//' Loglikelihood for a sample 
//' @param h positive parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' @param scale1 scale for observations
//' @param scale2 scale for censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMlognormal_UNIs(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  
  if(beta < 0) return NA_REAL;
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat s = (S * (-1)) * e;
  
  arma::mat mExp(size(S));
  
  double a = matrixMaxDiagonal_arma(S * (-1));
  
  int N{findN_arma(h, 1)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma(theVector, S, a, N);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{scale1[k] * pow(std::log(obs[k] + 1), beta)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * s;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta) + (beta - 1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{scale2[k] * pow(std::log(rcens[k] + 1), beta)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-loglogistic using Uniformization
//' 
//' Loglikelihood for a sample 
//' @param h positive parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' @param scale1 scale for observations
//' @param scale2 scale for censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMloglogistic_UNIs(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  
  if(beta[0] < 0 || beta[1] < 0) return NA_REAL;
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat s = (S * (-1)) * e;
  
  arma::mat mExp(size(S));
  
  double a = matrixMaxDiagonal_arma(S * (-1));
  
  int N{findN_arma(h, 1)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma(theVector, S, a, N);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{scale1[k] * std::log(pow(obs[k] / beta[0], beta[1]) + 1)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * s;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta[1]) - std::log(beta[0]) + (beta[1] - 1) * (std::log(obs[k]) - std::log(beta[0])) - std::log(pow(obs[k] / beta[0], beta[1]) + 1));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{scale2[k] * std::log(pow(rcens[k] / beta[0], beta[1]) + 1)};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}



//' Loglikelihood of matrix-Gompertz using Uniformization
//' 
//' Loglikelihood for a sample 
//' @param h positive parameter
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' @param scale1 scale for observations
//' @param scale2 scale for censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMgompertz_UNIs(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  
  if(beta < 0) return NA_REAL;
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat s = (S * (-1)) * e;
  
  arma::mat mExp(size(S));
  
  double a = matrixMaxDiagonal_arma(S * (-1));
  
  int N{findN_arma(h, 1)};
  
  std::vector<arma::mat> theVector;
  
  vectorOfMatrices_arma(theVector, S, a, N);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{scale1[k] * (exp(obs[k] * beta) - 1) / beta};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * s;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + obs[k] * beta);
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{scale2[k] * (exp(rcens[k] * beta) - 1) / beta};
    
    if (x * a <= 1.0) {
      mExp = matrixExpSum_arma(x, N, theVector, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      mExp = matrixExpSum_arma(x / pow(2.0, n), N, theVector, a);
      
      pow2Matrix_arma(n, mExp);
    }
    aux_mat = alpha.t() * mExp * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}







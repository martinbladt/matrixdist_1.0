#include <Rcpp.h>
using namespace Rcpp;
#include "matrix_functions.h"
#include "exp_arm.h"
#include "distributions.h"


//' Default size of the steps in the RK
//' 
//' Computes the default step length for a matrix \code{S} to be employed in the RK method
//' @param S sub-intensity matrix
//' @return The step length for \code{S}
//' 
// [[Rcpp::export]]
double default_step_length(const NumericMatrix & S) {
  double h{-0.1 / S(0,0)};
  
  for (int i{1}; i < S.nrow(); ++i) {
    if (h > -0.1 / S(i,i)) {
      h = -0.1 / S(i,i);
    }
  }
  return h;
}


//' Runge Kutta for the calculation of the a,b and c vectors in a EM step
//' 
//' Performce the RK of forth order
//' @param avector the a vector
//' @param bvector the b vector 
//' @param cmatrix the c matrix
//' @param dt the increment
//' @param h step-length
//' @param S sub-intensity
//' @param t exit rates 
//' 
// [[Rcpp::export]]
void runge_kutta(NumericMatrix & avector, NumericMatrix & bvector, NumericMatrix & cmatrix, double dt, double h, const NumericMatrix & S, const NumericMatrix & t) {
  int p{S.nrow()};
  int j{};
  int m{};
  double eps{};
  double sum{};
  
  int i{};
  i = dt / h;
  
  double h2{};
  h2 = dt / (i + 1);
  
  NumericMatrix ka(4,p);
  NumericMatrix kb(4,p);
  NumericMatrix kc1(p,p);
  NumericMatrix kc2(p,p);
  NumericMatrix kc3(p,p);
  NumericMatrix kc4(p,p);
  
  for (eps = 0; eps <= dt - h2 / 2; eps += h2) {
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(j,i) * avector(0,j);
      }
      ka(0,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(j,i) * (avector(0,j) + ka(0,j) / 2);
      }
      ka(1,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(j,i) * (avector(0,j) + ka(1,j) / 2);
      }
      ka(2,i) = h2 * sum;
    }
    for (i=0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j)
      {
        sum += S(j,i) * (avector(0,j) + ka(2,j));
      }
      ka(3,i) = h2 * sum;
    }
    
    
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(i,j) * bvector(j,0);
      }
      kb(0,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(i,j) * (bvector(j,0) + kb(0,j) / 2);
      }
      kb(1,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(i,j) * (bvector(j,0) + kb(1,j) / 2);
      }
      kb(2,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(i,j) * (bvector(j,0) + kb(2,j));
      }
      kb(3,i) = h2 * sum;
    }
    
    for (m = 0; m < p; ++m) {
      for (i = 0; i < p; ++i) {
        sum = t(m,0) * avector(0,i);
        for (j = 0; j < p; ++j) {
          sum += S(m,j) * cmatrix(j,i);
        }
        kc1(m,i) = h2 * sum;
      }
    }
    for (m = 0; m < p; ++m) {
      for (i = 0; i < p; ++i) {
        sum = t(m,0) * (avector(0,i) + ka(0,i) / 2);
        for (j = 0; j < p; ++j) {
          sum += S(m,j) * (cmatrix(j,i) + kc1(j,i) / 2);
        }
        kc2(m,i) = h2 * sum;
      }
    }
    for (m = 0; m < p; ++m) {
      for (i = 0; i < p; ++i) {
        sum = t(m,0) * (avector(0,i) + ka(1,i) / 2);
        for (j = 0; j < p; ++j) {
          sum += S(m,j) * (cmatrix(j,i) + kc2(j,i) / 2);
        }
        kc3(m,i) = h2 * sum;
      }
    }
    for (m = 0; m < p; ++m) {
      for (i = 0; i < p; ++i) {
        sum = t(m,0) * (avector(0,i) + ka(2,i));
        for (j = 0; j < p; ++j) {
          sum += S(m,j) * (cmatrix(j,i) + kc3(j,i));
        }
        kc4(m,i) = h2 * sum;
      }
    }
    
    for (i = 0; i < p; ++i) {
      avector(0,i) += (ka(0,i) + 2 * ka(1,i) + 2 * ka(2,i) + ka(3,i)) / 6;
      bvector(i,0) += (kb(0,i) + 2 * kb(1,i) + 2 * kb(2,i) + kb(3,i)) / 6;
      for (j = 0; j < p; ++j) {
        cmatrix(i,j) +=(kc1(i,j) + 2 * kc2(i,j) + 2 * kc3(i,j) + kc4(i,j)) / 6;
      }
    }
  }
}



//' EM step using Runge Kutta
//' 
//' Computes one step of the EM algorithm by using a Runge-Kutta method of 4th order
//' @param h step-length
//' @param alpha initial probalities
//' @param S sub-intensity
//' @param obs the observations
//' @param weight the weights for the observations
//' @param rcens censored observations
//' @param rcweight the weights for the censored observations
//' 
// [[Rcpp::export]]
void EMstep_RK(double h, NumericVector & alpha, NumericMatrix & S, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  long p{S.nrow()};
  
  NumericMatrix m_alpha(1, p, alpha.begin()); //Matrix version of alpha for computations
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  NumericMatrix Bmean(p,1);
  NumericMatrix Zmean(p,1);
  NumericMatrix Nmean(p,p + 1);
  
  NumericMatrix avector(1,p); 
  NumericMatrix bvector(p,1);
  NumericMatrix cmatrix(p,p);
  
  // initial conditions
  avector = clone(m_alpha);
  bvector = clone(t);
  
  double dt{0.0};
  if (obs.size()>0) {
    dt = obs[0];
  }
  
  double SumOfWeights{0.0};
  double density{0.0};
  
  // E step
  //  Uncensored data
  for (int k{0}; k < obs.size(); ++k) {
    
    SumOfWeights += weight[k];
    
    runge_kutta(avector, bvector, cmatrix, dt, h, S, t);
    density = matrix_product(m_alpha, bvector)(0,0);
    
    // E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha[i] * bvector(i,0) * weight[k] / density;
      Nmean(i,p) += avector(0,i) * t(i,0) * weight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * weight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += S(i,j) * cmatrix(j,i) * weight[k] / density;
      }
    }
    if (k < obs.size() - 1) {
      dt = obs[k + 1] - obs[k];
    }
  }
  
  //  Right-Censored Data
  double SumOfCensored{0.0};
  if (rcens.size() > 0) {
    dt = rcens[0];
    cmatrix = clone(cmatrix * 0);
    avector = clone(m_alpha);
    bvector = clone(e);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    
    SumOfCensored += rcweight[k];
    
    runge_kutta(avector, bvector, cmatrix, dt, h, S, e);
    density = matrix_product(m_alpha, bvector)(0,0);
    
    //E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha[i] * bvector(i,0) * rcweight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * rcweight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += S(i,j) * cmatrix(j,i) * rcweight[k] / density;
      }
    }
    if (k < rcens.size() - 1) {
      dt = rcens[k + 1] - rcens[k];
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


//' Runge Kutta for the calculation of the a vectors in a EM step 
//' 
//' Can be used for the loglikelihood
//' @param avector the a vector
//' @param dt increment
//' @param h step-length
//' @param S sub-intensity
//' 
// [[Rcpp::export]]
void a_rungekutta(NumericMatrix & avector, double dt, double h, const NumericMatrix & S) {
  long p{S.nrow()};
  int j{};
  double eps{};
  double sum{};
  
  int i{};
  i = dt / h;
  
  double h2{};
  h2 = dt / (i + 1);
  
  NumericMatrix ka(4,p);
  
  for (eps = 0; eps <= dt - h2 / 2; eps += h2) {
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(j,i) * avector(0,j);
      }
      ka(0,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(j,i) * (avector(0,j) + ka(0,j) / 2);
      }
      ka(1,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(j,i) * (avector(0,j) + ka(1,j) / 2);
      }
      ka(2,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(j,i) * (avector(0,j) + ka(2,j));
      }
      ka(3,i) = h2 * sum;
    }
    
    for (i = 0; i < p; ++i) {
      avector(0,i) += (ka(0,i) + 2 * ka(1,i) + 2 * ka(2,i) + ka(3,i)) / 6;
    }
  }
}


//' Loglikelihood using RK
//' 
//' Loglikelihood for a sample 
//' @param h step-length
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodPH_RK(double h, NumericVector & alpha, NumericMatrix & S, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  long p{S.nrow()};
  NumericMatrix m_alpha(1,p, alpha.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_alpha);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = obs[0];
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * std::log(density);
    if (k < obs.size() - 1) {dt = obs[k + 1] - obs[k]; }
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = rcens[0];
    avector = clone(m_alpha);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1) {dt = rcens[k + 1] - rcens[k];}
  }
  
  return logLh;
}


//' Loglikelihood of matrix Weibull using RK
//' 
//' Loglikelihood for a sample 
//' @param h step-length
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMweibull_RK(double h, NumericVector & alpha, NumericMatrix & S, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  
  if(beta < 0) return NA_REAL;
  
  long p{S.nrow()};
  NumericMatrix m_alpha(1,p, alpha.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_alpha);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = pow(obs[0], beta);
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta) + (beta -1) * std::log(obs[k]));
    if (k < obs.size() - 1) {dt = pow(obs[k + 1], beta) - pow(obs[k], beta);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = pow(rcens[0], beta);
    avector = clone(m_alpha);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1) {dt = pow(rcens[k + 1], beta) - pow(rcens[k], beta);}
  }
  
  return logLh;
}


//' Loglikelihood of matrix Pareto using RK
//' 
//' Loglikelihood for a sample 
//' @param h step-length
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMpareto_RK(double h, NumericVector & alpha, NumericMatrix & S, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  
  if(beta < 0) return NA_REAL;
  
  long p{S.nrow()};
  NumericMatrix m_alpha(1,p, alpha.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_alpha);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = std::log(obs[0] / beta + 1);
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (std::log(density) - std::log(obs[k] + beta));
    if (k < obs.size() - 1) {dt = std::log(obs[k + 1] / beta + 1) - std::log(obs[k] / beta + 1);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = std::log(rcens[0] / beta + 1);
    avector = clone(m_alpha);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1) {dt = std::log(rcens[k + 1] / beta + 1) - std::log(rcens[k] / beta + 1);}
  }
  
  return logLh;
}

//' Loglikelihood of matrix LogNormal using RK
//' 
//' Loglikelihood for a sample 
//' @param h step-length
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMlognormal_RK(double h, NumericVector & alpha, NumericMatrix & S, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  
  if(beta < 0) return NA_REAL;
  
  long p{S.nrow()};
  NumericMatrix m_alpha(1,p, alpha.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_alpha);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};

  // Non censored data
  if (obs.size() > 0) {
    dt = pow(std::log(obs[0] + 1), beta);
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta) + (beta -1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
    if (k < obs.size() - 1) {dt = pow(std::log(obs[k + 1] + 1), beta) - pow(std::log(obs[k] + 1), beta);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = pow(std::log(rcens[0] + 1), beta);
    avector = clone(m_alpha);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1) {dt = pow(std::log(rcens[k + 1] + 1), beta) - pow(std::log(rcens[k] + 1), beta);}
  }
  
  return logLh;
}

//' Loglikelihood of matrix Log-Logistic using RK
//' 
//' Loglikelihood for a sample 
//' @param h step-length
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMloglogistic_RK(double h, NumericVector & alpha, NumericMatrix & S, NumericVector beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  
  if(beta[0] < 0 || beta[1] < 0) return NA_REAL;
  
  
  long p{S.nrow()};
  NumericMatrix m_alpha(1,p, alpha.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_alpha);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = std::log(pow(obs[0] / beta[0], beta[1]) + 1);
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta[1]) - std::log(beta[0]) + (beta[1] - 1) * (std::log(obs[k]) - std::log(beta[0])) - std::log(pow(obs[k] / beta[0], beta[1]) + 1));
    if (k < obs.size() - 1) {dt = std::log(pow(obs[k + 1] / beta[0], beta[1]) + 1) - std::log(pow(obs[k] / beta[0], beta[1]) + 1);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = std::log(pow(rcens[0] / beta[0], beta[1]) + 1);
    avector = clone(m_alpha);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1) {dt = std::log(pow(rcens[k + 1] / beta[0], beta[1]) + 1) - std::log(pow(rcens[k] / beta[0], beta[1]) + 1);}
  }
  
  return logLh;
}

//' Loglikelihood of matrix Gompertz using RK
//' 
//' Loglikelihood for a sample 
//' @param h step-length
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMgompertz_RK(double h, NumericVector & alpha, NumericMatrix & S, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  
  if(beta < 0) return NA_REAL;
  
  long p{S.nrow()};
  NumericMatrix m_alpha(1,p, alpha.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_alpha);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = (exp(obs[0] * beta) - 1) / beta;
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (std::log(density) + obs[k] * beta);
    if (k < obs.size() - 1) {dt = (exp(obs[k + 1] * beta) - 1) / beta - (exp(obs[k] * beta) - 1) / beta;}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = (exp(rcens[0] * beta) - 1) / beta;
    avector = clone(m_alpha);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1) {dt = (exp(rcens[k + 1] * beta) - 1) / beta - (exp(rcens[k] * beta) - 1) / beta;}
  }
  
  return logLh;
}


//' Loglikelihood of matrix GEV using RK
//' 
//' Loglikelihood for a sample 
//' @param h step-length
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double logLikelihoodMgev_RK(double h, NumericVector & alpha, NumericMatrix & S, NumericVector beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  
  if(beta[1] < 0) return NA_REAL;
  
  long p{S.nrow()};
  NumericMatrix m_alpha(1,p, alpha.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_alpha);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  int N = static_cast<int>(obs.size());
  
  if (beta[2] == 0) {
    // Non censored data
    if (N > 0) {
      dt = exp(-(obs[N - 1] - beta[0]) / beta[1]);
    }
    for (int k{1}; k <= N; ++k) {
      a_rungekutta(avector, dt, h, S);
      density = matrix_product(avector, t)(0,0);
      logLh += weight[N - k] * (std::log(density) - std::log(beta[1]) - (obs[N - k] - beta[0]) / beta[1]);
      if (k < N) {dt = exp(-(obs[N - k - 1] - beta[0]) / beta[1]) - exp(-(obs[N - k] - beta[0]) / beta[1]);}
    }
    //Right censored data
    N = rcens.size();
    if (N > 0) {
      dt = exp(-(rcens[N - 1] - beta[0]) / beta[1]);
      avector = clone(m_alpha);
    }
    for (int k{1}; k <= N; ++k) {
      a_rungekutta(avector, dt, h, S);
      density = matrix_product(avector, e)(0,0);
      logLh += rcweight[N - k] * std::log(density);
      if (k < N) { dt = exp(-(rcens[N - k - 1] - beta[0]) / beta[1]) - exp(-(rcens[N - k] - beta[0]) / beta[1]);}
    }
  }
  else {
    // Non censored data
    if (N > 0) {
      dt = pow(1 + (beta[2] / beta[1]) * (obs[N - 1] - beta[0]) , - 1 / beta[2]);
    }
    for (int k{1}; k <= N; ++k) {
      a_rungekutta(avector, dt, h, S);
      density = matrix_product(avector, t)(0,0);
      logLh += weight[N - k] * (std::log(density) - std::log(beta[1]) - (1 + 1 / beta[2]) * std::log(1 + (beta[2] / beta[1]) * (obs[N - k] - beta[0])));
      if (k < N) {dt = pow(1 + (beta[2] / beta[1]) * (obs[N - k - 1] - beta[0]) , - 1 / beta[2]) - pow(1 + (beta[2] / beta[1]) * (obs[N - k] - beta[0]) , - 1 / beta[2]);}
      
    }
    //Right censored data
    N = rcens.size();
    if (N > 0) {
      dt = pow(1 + (beta[2] / beta[1]) * (rcens[N - 1] - beta[0]) , - 1 / beta[2]);
      avector = clone(m_alpha);
    }
    for (int k{1}; k <= N; ++k) {
      a_rungekutta(avector, dt, h, S);
      density = matrix_product(avector, e)(0,0);
      logLh += rcweight[N - k] * std::log(density);
      if (k < N) {dt = pow(1 + (beta[2] / beta[1]) * (rcens[N - k - 1] - beta[0]) , - 1 / beta[2]) - pow(1 + (beta[2] / beta[1]) * (rcens[N - k] - beta[0]) , - 1 / beta[2]);}
    }
  }
  return logLh;
}


//' Applies the inverse of the GEV but giving back the vector in reverse order
//' 
//' Used for EM step
//' @param observations the observations
//' @param weights weithgs of the observations
//' @param beta parameters of the GEV
//' 
// [[Rcpp::export]]
List reversTransformData(const NumericVector & observations, const NumericVector & weights, const NumericVector & beta) {
  int N = static_cast<int>(observations.size());
  NumericVector SransformObs(N);
  NumericVector SransWeights(N);
  if (beta[2] == 0) { // Gumbel
    for (int i{0}; i < N; ++i) {
      SransformObs[i] = exp( -(observations[N - i - 1] - beta[0]) / beta[1]) ;
      SransWeights[i] = weights[N - i - 1];
    }
  }
  else { // GEVD
    for (int i{0}; i < N; ++i) {
      SransformObs[i] = pow( 1 + beta[2] * (observations[N - i - 1] - beta[0]) / beta[1] , -1 / beta[2]);
      SransWeights[i] = weights[N - i - 1];
    }
  }
  
  List L = List::create(Named("obs") = SransformObs, _["weight"] = SransWeights);
  
  return L;
}

//' Derivative of matrix Weibull
//' 
//' Can be used to increase performance
//' @param h step-length
//' @param alpha initial probabilities
//' @param S sub-intensity
//' @param beta parameter of transformation
//' @param obs the observations
//' @param weight weight of the observations
//' @param rcens censored observations
//' @param rcweight weight of the censored observations
//' 
// [[Rcpp::export]]
double derivativeMatrixweibull(double h, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, NumericVector & alpha, NumericMatrix & S,  double beta) {
  long p{S.nrow()};
  NumericMatrix m_alpha(1,p, alpha.begin());
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  NumericVector aux_vet(1);
  
  double logLh{0.0};
  for (int k{0}; k < obs.size(); ++k) {
    aux_vet[0] = pow(obs[k], beta);
    logLh +=  weight[k] * ( (matrix_product(m_alpha, matrix_product(matrix_exponential(S * pow(obs[k], beta)), matrix_product(S, t)))(0,0) * std::log(obs[k]) * pow(obs[k], beta) ) / mweibullden(aux_vet, alpha, S, beta)[0] + 1 / beta + std::log(obs[k]) );
  }
  return logLh;
  
  for (int k{0}; k < rcens.size(); ++k) {
    aux_vet[0] = pow(rcens[k], beta);
    logLh +=  rcweight[k] * ( (matrix_product(m_alpha, matrix_product(matrix_exponential(S * pow(obs[k], beta)), matrix_product(S, e)))(0,0) * std::log(obs[k]) * pow(obs[k], beta) ) / mweibullcdf(aux_vet, alpha, S, beta, false)[0] );
  }
  return logLh;
  
}


////////////////////////////////////////////
// Scaled versions of loglikelihoods (for regression):
////////////////////////////////////////////


//' Loglikelihood of PH using RK
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
double logLikelihoodPH_RKs(double h, NumericVector & alpha, NumericMatrix & S, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{S.nrow()};
  NumericMatrix m_alpha(1,p, alpha.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_alpha);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * obs[0];
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]));
    if (k < obs.size() - 1){dt = scale1[k + 1] * obs[k + 1] - scale1[k] * obs[k];}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * rcens[0];
    avector = clone(m_alpha);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * rcens[k + 1] - scale2[k] * rcens[k];}
  }
  
  return logLh;
}


//' Loglikelihood of matrix-Weibull using RK
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
double logLikelihoodMweibull_RKs(double h, NumericVector & alpha, NumericMatrix & S, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{S.nrow()};
  NumericMatrix m_alpha(1,p, alpha.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_alpha);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * pow(obs[0], beta);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta) + (beta -1) * std::log(obs[k]));
    if (k < obs.size() - 1){dt = scale1[k + 1] * pow(obs[k + 1], beta) - scale1[k] * pow(obs[k], beta);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * pow(rcens[0], beta);
    avector = clone(m_alpha);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * pow(rcens[k + 1], beta) - scale2[k] * pow(rcens[k], beta);}
  }
  
  return logLh;
}


//' Loglikelihood of matrix-Pareto using RK
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
double logLikelihoodMpareto_RKs(double h, NumericVector & alpha, NumericMatrix & S, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{S.nrow()};
  NumericMatrix m_alpha(1,p, alpha.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_alpha);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * std::log(obs[0] / beta + 1);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) - std::log(obs[k] + beta));
    if (k < obs.size() - 1){dt = scale1[k + 1] * std::log(obs[k + 1] / beta + 1) - scale1[k] * std::log(obs[k] / beta + 1);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * std::log(rcens[0] / beta + 1);
    avector = clone(m_alpha);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * std::log(rcens[k + 1] / beta + 1) - scale2[k] * std::log(rcens[k] / beta + 1);}
  }
  
  return logLh;
}


//' Loglikelihood of matrix-lognormal using RK
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
double logLikelihoodMlognormal_RKs(double h, NumericVector & alpha, NumericMatrix & S, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{S.nrow()};
  NumericMatrix m_alpha(1,p, alpha.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_alpha);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * pow(std::log(obs[0] + 1), beta);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta) + (beta -1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
    if (k < obs.size() - 1){dt = scale1[k + 1] * pow(std::log(obs[k + 1] + 1), beta) - scale1[k] * pow(std::log(obs[k] + 1), beta);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * pow(std::log(rcens[0] + 1), beta);
    avector = clone(m_alpha);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * pow(std::log(rcens[k + 1] + 1), beta) - scale2[k] * pow(std::log(rcens[k] + 1), beta);}
  }
  
  return logLh;
}



//' Loglikelihood of matrix-loglogistic using RK
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
double logLikelihoodMloglogistic_RKs(double h, NumericVector & alpha, NumericMatrix & S, NumericVector beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{S.nrow()};
  NumericMatrix m_alpha(1,p, alpha.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_alpha);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * std::log(pow(obs[0] / beta[0], beta[1]) + 1);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta[1]) - std::log(beta[0]) + (beta[1] - 1) * (std::log(obs[k]) - std::log(beta[0])) - std::log(pow(obs[k] / beta[0], beta[1]) + 1));
    if (k < obs.size() - 1){dt = scale1[k + 1] * std::log(pow(obs[k + 1] / beta[0], beta[1]) + 1) - scale1[k] * std::log(pow(obs[k] / beta[0], beta[1]) + 1);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * std::log(pow(rcens[0] / beta[0], beta[1]) + 1);
    avector = clone(m_alpha);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * std::log(pow(rcens[k + 1] / beta[0], beta[1]) + 1) - scale2[k] * std::log(pow(rcens[k] / beta[0], beta[1]) + 1);}
  }
  
  return logLh;
}



//' Loglikelihood of matrix-Gompertz using RK
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
double logLikelihoodMgompertz_RKs(double h, NumericVector & alpha, NumericMatrix & S, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{S.nrow()};
  NumericMatrix m_alpha(1,p, alpha.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(S * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_alpha);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * (exp(obs[0] * beta) - 1) / beta;
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + obs[k] * beta);
    if (k < obs.size() - 1){dt = scale1[k + 1] * (exp(obs[k + 1] * beta) - 1) / beta - scale1[k] * (exp(obs[k] * beta) - 1) / beta;}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * (exp(rcens[0] * beta) - 1) / beta;
    avector = clone(m_alpha);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * (exp(rcens[k + 1] * beta) - 1) / beta - scale2[k] * (exp(rcens[k] * beta) - 1) / beta;}
  }
  
  return logLh;
}





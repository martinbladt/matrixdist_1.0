#include <Rcpp.h>
using namespace Rcpp;
#include "matrix_functions.h"


//' Default size of the steps in the RK
// [[Rcpp::export]]
double default_step_length(const NumericMatrix & T) {
  double h{-0.1 / T(0,0)};
  
  for (int i{1}; i < T.nrow(); ++i) {
    if (h > -0.1 / T(i,i)) {
      h = -0.1 / T(i,i);
    }
  }
  return h;
}


//' Runge Kutta for the calculation of the a,b and c vectors in a EM step
//' 
//' I may need to change the type of avector and bvector, depending on how I call them in the EM step
// [[Rcpp::export]]
void runge_kutta(NumericMatrix & avector, NumericMatrix & bvector, NumericMatrix & cmatrix, double dt, double h, const NumericMatrix & T, const NumericMatrix & t) {
  int p{T.nrow()};
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
        sum += T(j,i) * avector(0,j);
      }
      ka(0,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += T(j,i) * (avector(0,j) + ka(0,j) / 2);
      }
      ka(1,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += T(j,i) * (avector(0,j) + ka(1,j) / 2);
      }
      ka(2,i) = h2 * sum;
    }
    for (i=0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j)
      {
        sum += T(j,i) * (avector(0,j) + ka(2,j));
      }
      ka(3,i) = h2 * sum;
    }
    
    
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += T(i,j) * bvector(j,0);
      }
      kb(0,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += T(i,j) * (bvector(j,0) + kb(0,j) / 2);
      }
      kb(1,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += T(i,j) * (bvector(j,0) + kb(1,j) / 2);
      }
      kb(2,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += T(i,j) * (bvector(j,0) + kb(2,j));
      }
      kb(3,i) = h2 * sum;
    }
    
    for (m = 0; m < p; ++m) {
      for (i = 0; i < p; ++i) {
        sum = t(m,0) * avector(0,i);
        for (j = 0; j < p; ++j) {
          sum += T(m,j) * cmatrix(j,i);
        }
        kc1(m,i) = h2 * sum;
      }
    }
    for (m = 0; m < p; ++m) {
      for (i = 0; i < p; ++i) {
        sum = t(m,0) * (avector(0,i) + ka(0,i) / 2);
        for (j = 0; j < p; ++j) {
          sum += T(m,j) * (cmatrix(j,i) + kc1(j,i) / 2);
        }
        kc2(m,i) = h2 * sum;
      }
    }
    for (m = 0; m < p; ++m) {
      for (i = 0; i < p; ++i) {
        sum = t(m,0) * (avector(0,i) + ka(1,i) / 2);
        for (j = 0; j < p; ++j) {
          sum += T(m,j) * (cmatrix(j,i) + kc2(j,i) / 2);
        }
        kc3(m,i) = h2 * sum;
      }
    }
    for (m = 0; m < p; ++m) {
      for (i = 0; i < p; ++i) {
        sum = t(m,0) * (avector(0,i) + ka(2,i));
        for (j = 0; j < p; ++j) {
          sum += T(m,j) * (cmatrix(j,i) + kc3(j,i));
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
// [[Rcpp::export]]
void EMstep_RK(double h, NumericVector & pi, NumericMatrix & T, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  long p{T.nrow()};
  
  NumericMatrix m_pi(1, p, pi.begin()); //Matrix version of pi for computations
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  NumericMatrix Bmean(p,1);
  NumericMatrix Zmean(p,1);
  NumericMatrix Nmean(p,p + 1);
  
  NumericMatrix avector(1,p); 
  NumericMatrix bvector(p,1);
  NumericMatrix cmatrix(p,p);
  
  // initial conditions
  avector = clone(m_pi);
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
    
    runge_kutta(avector, bvector, cmatrix, dt, h, T, t);
    density = matrix_product(m_pi, bvector)(0,0);
    
    // E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += pi[i] * bvector(i,0) * weight[k] / density;
      Nmean(i,p) += avector(0,i) * t(i,0) * weight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * weight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += T(i,j) * cmatrix(j,i) * weight[k] / density;
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
    avector = clone(m_pi);
    bvector = clone(e);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    
    SumOfCensored += rcweight[k];
    
    runge_kutta(avector, bvector, cmatrix, dt, h, T, e);
    density = matrix_product(m_pi, bvector)(0,0);
    
    //E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += pi[i] * bvector(i,0) * rcweight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * rcweight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += T(i,j) * cmatrix(j,i) * rcweight[k] / density;
      }
    }
    if (k < rcens.size() - 1) {
      dt = rcens[k + 1] - rcens[k];
    }
  }
  
  // M step
  for (int i{0}; i < p; ++i) {
    pi[i] = Bmean(i,0) / (SumOfWeights + SumOfCensored);
    if (pi[i] < 0) {
      pi[i] = 0;
    }
    t(i,0) = Nmean(i,p) / Zmean(i,0);
    if (t(i,0) < 0) {
      t(i,0) = 0;
    }
    T(i,i) = -t(i,0);
    for (int j{0}; j < p; ++j) {
      if (i != j) {
        T(i,j) = Nmean(i,j) / Zmean(i,0);
        if (T(i,j) < 0) {
          T(i,j) = 0;
        }
        T(i,i) -= T(i,j);
      }
    }
  }
}


//' Runge Kutta for the calculation of the a vectors in a EM step - Can be used for the loglikelihood
//' 
//' I may need to change the type of avector
// [[Rcpp::export]]
void a_rungekutta(NumericMatrix & avector, double dt, double h, const NumericMatrix & T) {
  long p{T.nrow()};
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
        sum += T(j,i) * avector(0,j);
      }
      ka(0,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += T(j,i) * (avector(0,j) + ka(0,j) / 2);
      }
      ka(1,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += T(j,i) * (avector(0,j) + ka(1,j) / 2);
      }
      ka(2,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += T(j,i) * (avector(0,j) + ka(2,j));
      }
      ka(3,i) = h2 * sum;
    }
    
    for (i = 0; i < p; ++i) {
      avector(0,i) += (ka(0,i) + 2 * ka(1,i) + 2 * ka(2,i) + ka(3,i)) / 6;
    }
  }
}


//' Loglikelihood using RK
// [[Rcpp::export]]
double logLikelihoodPH_RK(double h, NumericVector & pi, NumericMatrix & T, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = obs[0];
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * log(density);
    if (k < obs.size() - 1) {dt = obs[k + 1] - obs[k]; }
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = rcens[0];
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1) {dt = rcens[k + 1] - rcens[k];}
  }
  
  return logLh;
}



//' Loglikelihood IPH using RK and g as an input
//' One needs to be careful with the GEV since it is decreasing
//' It is slower than using the density directly - Perhaps is the iteration with R 
// [[Rcpp::export]]
double logLikelihoodIPH_RK(double h, NumericVector & pi, NumericMatrix & T, Function g, Function g_inv, Function lambda, NumericVector beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  NumericVector g_inv_val;
  NumericVector g_inv_val2;
  NumericVector lambda_val;
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    g_inv_val = g_inv(obs[0], beta);
    dt = g_inv_val[0];
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    lambda_val = lambda(obs[k], beta);
    logLh += weight[k] * (log(density) + log(lambda_val[0]));
    g_inv_val = g_inv(obs[k], beta);
    if (k < obs.size() - 1) {
      g_inv_val2 = g_inv(obs[k + 1], beta);
      dt = g_inv_val2[0] - g_inv_val[0];}
  }
  //Right censored data
  if (rcens.size() > 0) {
    g_inv_val = g_inv(rcens[0], beta);
    dt = g_inv_val[0];
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    g_inv_val = g_inv(rcens[k], beta);
    if (k < rcens.size() - 1) {
      g_inv_val2 = g_inv(rcens[k + 1], beta);
      dt = g_inv_val2[0] - g_inv_val[0];
    }
  }
  
  return logLh;
}



//' Loglikelihood of matrix Weibull using RK
//' This is the fastest option
// [[Rcpp::export]]
double logLikelihoodMWeib_RK(double h, NumericVector & pi, NumericMatrix & T, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = pow(obs[0], beta);
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) + log(beta) + (beta -1) * log(obs[k]));
    if (k < obs.size() - 1) {dt = pow(obs[k + 1], beta) - pow(obs[k], beta);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = pow(rcens[0], beta);
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1) {dt = pow(rcens[k + 1], beta) - pow(rcens[k], beta);}
  }
  
  return logLh;
}



//' Loglikelihood of matrix Pareto using RK
// [[Rcpp::export]]
double logLikelihoodMPar_RK(double h, NumericVector & pi, NumericMatrix & T, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = log(obs[0] / beta + 1);
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) - log(obs[k] + beta));
    if (k < obs.size() - 1) {dt = log(obs[k + 1] / beta + 1) - log(obs[k] / beta + 1);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = log(rcens[0] / beta + 1);
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1) {dt = log(rcens[k + 1] / beta + 1) - log(rcens[k] / beta + 1);}
  }
  
  return logLh;
}


//' Loglikelihood of matrix Log-Logistic using RK
// [[Rcpp::export]]
double logLikelihoodMLogLogistic_RK(double h, NumericVector & pi, NumericMatrix & T, NumericVector beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = log(pow(obs[0] / beta[0], beta[1]) + 1);
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) + log(beta[1]) - log(beta[0]) + (beta[1] - 1) * (log(obs[k]) - log(beta[0])) - log(pow(obs[k] / beta[0], beta[1]) + 1));
    if (k < obs.size() - 1) {dt = log(pow(obs[k + 1] / beta[0], beta[1]) + 1) - log(pow(obs[k] / beta[0], beta[1]) + 1);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = log(pow(rcens[0] / beta[0], beta[1]) + 1);
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1) {dt = log(pow(rcens[k + 1] / beta[0], beta[1]) + 1) - log(pow(rcens[k] / beta[0], beta[1]) + 1);}
  }
  
  return logLh;
}

//' Loglikelihood of matrix Gompertz using RK
// [[Rcpp::export]]
double logLikelihoodMGomp_RK(double h, NumericVector & pi, NumericMatrix & T, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = (exp(obs[0] * beta) - 1) / beta;
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) + obs[k] * beta);
    if (k < obs.size() - 1) {dt = (exp(obs[k + 1] * beta) - 1) / beta - (exp(obs[k] * beta) - 1) / beta;}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = (exp(rcens[0] * beta) - 1) / beta;
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1) {dt = (exp(rcens[k + 1] * beta) - 1) / beta - (exp(rcens[k] * beta) - 1) / beta;}
  }
  
  return logLh;
}




//' Loglikelihood of matrix GEV using RK
//' I am assuming that the sample is given in an increasing order
// [[Rcpp::export]]
double logLikelihoodMGEV_RK(double h, NumericVector & pi, NumericMatrix & T, NumericVector beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  long N{obs.size()};
  
  if (beta[2] == 0) {
    // Non censored data
    if (N > 0) {
      dt = exp(-(obs[N - 1] - beta[0]) / beta[1]);
    }
    for (int k{1}; k <= N; ++k) {
      a_rungekutta(avector, dt, h, T);
      density = matrix_product(avector, t)(0,0);
      logLh += weight[N - k] * (log(density) - log(beta[1]) - (obs[N - k] - beta[0]) / beta[1]);
      if (k < N) {dt = exp(-(obs[N - k - 1] - beta[0]) / beta[1]) - exp(-(obs[N - k] - beta[0]) / beta[1]);}
    }
    //Right censored data
    N = rcens.size();
    if (N > 0) {
      dt = exp(-(rcens[N - 1] - beta[0]) / beta[1]);
      avector = clone(m_pi);
    }
    for (int k{1}; k <= N; ++k) {
      a_rungekutta(avector, dt, h, T);
      density = matrix_product(avector, e)(0,0);
      logLh += rcweight[N - k] * log(density);
      if (k < N) { dt = exp(-(rcens[N - k - 1] - beta[0]) / beta[1]) - exp(-(rcens[N - k] - beta[0]) / beta[1]);}
    }
  }
  else {
    // Non censored data
    if (N > 0) {
      dt = pow(1 + (beta[2] / beta[1]) * (obs[N - 1] - beta[0]) , - 1 / beta[2]);
    }
    for (int k{1}; k <= N; ++k) {
      a_rungekutta(avector, dt, h, T);
      density = matrix_product(avector, t)(0,0);
      logLh += weight[N - k] * (log(density) - log(beta[1]) - (1 + 1 / beta[2]) * log(1 + (beta[2] / beta[1]) * (obs[N - k] - beta[0])));
      if (k < N) {dt = pow(1 + (beta[2] / beta[1]) * (obs[N - k - 1] - beta[0]) , - 1 / beta[2]) - pow(1 + (beta[2] / beta[1]) * (obs[N - k] - beta[0]) , - 1 / beta[2]);}
      
    }
    //Right censored data
    N = rcens.size();
    if (N > 0) {
      dt = pow(1 + (beta[2] / beta[1]) * (rcens[N - 1] - beta[0]) , - 1 / beta[2]);
      avector = clone(m_pi);
    }
    for (int k{1}; k <= N; ++k) {
      a_rungekutta(avector, dt, h, T);
      density = matrix_product(avector, e)(0,0);
      logLh += rcweight[N - k] * log(density);
      if (k < N) {dt = pow(1 + (beta[2] / beta[1]) * (rcens[N - k - 1] - beta[0]) , - 1 / beta[2]) - pow(1 + (beta[2] / beta[1]) * (rcens[N - k] - beta[0]) , - 1 / beta[2]);}
    }
  }
  return logLh;
}


//' Applies the inverse of the GEV but giving back the vector in reverse order
// [[Rcpp::export]]
List reversTransformData(const NumericVector & observations, const NumericVector & weights, const NumericVector & beta) {
  long N{observations.size()};
  NumericVector TransformObs(N);
  NumericVector TransWeights(N);
  if (beta[2] == 0) { // Gumbel
    for (int i{0}; i < N; ++i) {
      TransformObs[i] = exp( -(observations[N - i - 1] - beta[0]) / beta[1]) ;
      TransWeights[i] = weights[N - i - 1];
    }
  }
  else { // GEVD
    for (int i{0}; i < N; ++i) {
      TransformObs[i] = pow( 1 + beta[2] * (observations[N - i - 1] - beta[0]) / beta[1] , -1 / beta[2]);
      TransWeights[i] = weights[N - i - 1];
    }
  }
  
  List L = List::create(Named("obs") = TransformObs, _["weight"] = TransWeights);
  
  return L;
}




//'  EM for a Bivariate PH fit
//'  
// [[Rcpp::export]]
void EMstep_bivph(const NumericMatrix & observations, const NumericVector & weights, NumericVector & alpha, NumericMatrix & T11, NumericMatrix & T12, NumericMatrix & T22) {
  long p1{T11.nrow()};
  long p2{T22.nrow()};
  long p{p1 + p2};
  
  double density{0};
  double aux;
  
  NumericMatrix m_alpha(1,p1, alpha.begin());
  
  NumericMatrix Bmean(p1,1);
  NumericMatrix Zmean(p,1);
  NumericMatrix Nmean(p,p + 1);
  
  NumericMatrix bmatrix1(p1,p1);
  NumericMatrix bmatrix2(p2,p2);
  NumericMatrix cmatrix1(p1,p1);
  NumericMatrix cmatrix2(p2,p2);
  NumericMatrix aux_exp1(p1,p1);
  NumericMatrix aux_exp2(p2,p2);
  
  NumericMatrix J1(2 * p1,2 * p1);
  NumericMatrix J2(2 * p2,2 * p2);
  
  
  NumericVector m_e(p2, 1);
  NumericMatrix e(p2, 1, m_e.begin());
  
  NumericMatrix exitvec = matrix_product(T22 * (-1), e);
  
  NumericMatrix auxMatrix1(p1,1);
  NumericMatrix auxMatrix2(p2,p1);
  NumericMatrix auxMatrix3(1,p2);
  
  double sumOfWeights{0.0};
  //E step
  for (int k{0}; k < observations.nrow(); ++k) {
    
    sumOfWeights += weights[k];
    
    aux_exp2 = matrix_exponential(T22 * observations(k,1));
    bmatrix1 = matrix_product(T12, matrix_product(aux_exp2, matrix_product(exitvec, m_alpha)));
    
    J1 = matrix_exponential(matrix_VanLoan(T11, T11, bmatrix1) * observations(k,0)); 
    
    for (int i{0}; i < p1; ++i) {
      for (int j{0}; j < p1; ++j) {
        aux_exp1(i,j) = J1(i,j);
        cmatrix1(i,j) = J1(i,j + p1);
      }
    }
    
    bmatrix2 = matrix_product(exitvec, matrix_product(m_alpha, matrix_product(aux_exp1, T12)));
    
    J2 = matrix_exponential(matrix_VanLoan(T22, T22, bmatrix2) * observations(k,1)); 
    
    for (int i{0}; i < p2; ++i) {
      for (int j{0}; j < p2; ++j) {
        cmatrix2(i,j) = J2(i,j + p2);
      }
    }
    density = matrix_product(m_alpha, matrix_product(aux_exp1, matrix_product(T12, matrix_product(aux_exp2, exitvec))))(0,0);
    
    
    //E-step
    auxMatrix1 = matrix_product(aux_exp1, matrix_product(T12, matrix_product(aux_exp2, exitvec)));
    auxMatrix2 =  matrix_product(aux_exp2, matrix_product(exitvec, matrix_product(m_alpha, aux_exp1)));
    for (int i{0}; i < p1; ++i) {
      aux = auxMatrix1(i,0);
      Bmean(i,0) += m_alpha(0,i) * aux * weights[k] / density;
      Zmean(i,0) += cmatrix1(i,i) * weights[k] / density;
      for (int j{0}; j < p1; ++j) {
        Nmean(i,j) += T11(i,j) * cmatrix1(j,i) * weights[k] / density;
      }
      for (int j{0}; j < p2; ++j) {
        aux = auxMatrix2(j,i);
        Nmean(i,j + p1) += T12(i,j) * aux * weights[k] / density;
      }
    }
    
    auxMatrix3 = matrix_product(m_alpha, matrix_product(aux_exp1, matrix_product(T12, aux_exp2)));
    for (int i{0}; i < p2; ++i) {
      Zmean(i + p1,0) += cmatrix2(i,i) * weights[k] / density;
      aux = auxMatrix3(0,i);
      Nmean(i + p1,p) += aux * exitvec(i,0) * weights[k] / density;
      for (int j{0}; j < p2; ++j){
        Nmean(i + p1,j + p1) += T22(i,j) * cmatrix2(j,i) * weights[k] / density;
      }
    }
  }
  
  // M step
  for (int i{0}; i < p1; ++i) {
    alpha[i] = Bmean(i,0) / sumOfWeights;
    if (alpha[i] < 0) {
      alpha[i] = 0;
    }
    T11(i,i) = 0;
    for (int j{0}; j < p1; ++j) {
      if (i != j) {
        T11(i,j) = Nmean(i,j) / Zmean(i,0);
        if (T11(i,j) < 0) {
          T11(i,j) = 0;
        }
        T11(i,i) -= T11(i,j);
      }
    }
    for (int j{0}; j < p2; ++j) {
      T12(i,j) = Nmean(i,j + p1) / Zmean(i,0);
      if (T12(i,j) < 0){
        T12(i,j) = 0;
      }
      T11(i,i) -= T12(i,j);
    }
  }
  for (int i{0}; i < p2; ++i) {
    exitvec(i,0) = Nmean(i + p1,p) / Zmean(i + p1,0);
    if (exitvec(i,0) < 0) {
      exitvec(i,0) = 0;
    }
    T22(i,i) = -exitvec(i,0);
    for (int j{0}; j < p2; ++j) {
      if (i != j) {
        T22(i,j) = Nmean(i + p1,j + p1) / Zmean(i + p1,0);
        if (T22(i,j) < 0){
          T22(i,j) = 0;
        }
        T22(i,i) -= T22(i,j);
      }
    }
  }
  
}



//' EM using Matlab algorithm for matrix exponential
//' 
//' This one is slower but dont requires to order the sample
// [[Rcpp::export]]
void EMstep(NumericVector & pi, NumericMatrix & T, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight) {
  long p{T.nrow()};
  
  NumericMatrix m_pi(1, p, pi.begin()); //Matrix version of pi for computations
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  NumericMatrix Bmean(p,1);
  NumericMatrix Zmean(p,1);
  NumericMatrix Nmean(p,p + 1);
  
  NumericMatrix avector(1,p);
  NumericMatrix bvector(p,1);
  NumericMatrix cmatrix(p,p);
  NumericMatrix aux_exp(p,p);
  
  NumericMatrix J(2 * p,2 * p);
  NumericMatrix tProductPi(p,p);
  tProductPi = matrix_product(t, m_pi);
  
  double SumOfWeights{0.0};
  double density{0.0};
  
  //E-step
  //  Unccensored data
  for (int k{0}; k < obs.size(); ++k) {
    SumOfWeights += weight[k];
    
    J = matrix_exponential(matrix_VanLoan(T, T, tProductPi) * obs[k]);
    
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }
    
    avector = matrix_product(m_pi, aux_exp);
    bvector = matrix_product(aux_exp, t);
    density = matrix_product(m_pi, bvector)(0,0);
    
    //E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += m_pi(0,i) * bvector(i,0) * weight[k] / density;
      Nmean(i,p) += avector(0,i) * t(i,0) * weight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * weight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += T(i,j) * cmatrix(j,i) * weight[k] / density;
      }
    }
  }
  //  Right-Censored Data
  double SumOfCensored{0.0};
  if (rcens.size() > 0) {
    tProductPi = matrix_product(e, m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    SumOfCensored += rcweight[k];
    
    J = matrix_exponential(matrix_VanLoan(T, T, tProductPi) * rcens[k]);
    
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }
    
    bvector = matrix_product(aux_exp, e);
    density = matrix_product(m_pi, bvector)(0,0);
    
    //E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += m_pi(0,i) * bvector(i,0) * rcweight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * rcweight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += T(i,j) * cmatrix(j,i) * rcweight[k] / density;
      }
    }
  }
  
  // M step
  for (int i{0}; i < p; ++i) {
    pi[i] = Bmean(i,0) / (SumOfWeights + SumOfCensored);
    if (pi[i] < 0) {
      pi[i] = 0;
    }
    t(i,0) = Nmean(i,p) / Zmean(i,0);
    if (t(i,0) < 0) {
      t(i,0) = 0;
    }
    T(i,i) = -t(i,0);
    for (int j{0}; j < p; ++j) {
      if (i != j) {
        T(i,j) = Nmean(i,j) / Zmean(i,0);
        if (T(i,j) < 0) {
          T(i,j) = 0;
        }
        T(i,i) -= T(i,j);
      }
    }
  }
}



//' Pi and T of a linear combination of a MPH*
//' 
//' @examples
//' pi <- c(0.15, 0.85, 0 ,0)
//' T11 <- matrix(c(c(-2,9),c(0,-11)), nrow = 2, ncol = 2)
//' T12 <- matrix(c(c(2,0),c(0,2)), nrow = 2, ncol = 2)
//' T22 <- matrix(c(c(-1,0),c(0.5,-5)), nrow = 2, ncol = 2)
//' T <- merge_matrices(T11, T12, T22)
//' R <- matrix(c(c(1,1,0,0), c(0,0,1,1)), ncol=2)
//' w1 <- c(1,0)
//' linear_combination(w1, pi, T, R)
//' w2 <- c(0,1)
//' linear_combination(w2, pi, T, R)
//' matrix(c(0.15, 0.85), ncol=2)%*%matrix_inverse(T11 * (-1))%*%T12
//' w3 <- c(1,1)
//' linear_combination(w3, pi, T, R)
// [[Rcpp::export]]
List linear_combination(NumericVector w, NumericVector pi, NumericMatrix T, NumericMatrix R) {
  long p{T.nrow()};
  
  NumericVector newStates;
  
  int NumZeros{0};
  IntegerVector deleteRows; //states to be deleted
  IntegerVector keepRows; //states to keep
  
  NumericMatrix m_w(w.size(), 1, w.begin());
  
  NumericMatrix Rw(matrix_product(R, m_w)); 
  
  NumericMatrix m_pi(1, pi.size(), pi.begin());
  
  for (int j{0}; j < p; ++j) {
    if (Rw(j,0) == 0) {
      deleteRows.push_back(j);
      ++NumZeros;
    }
    else {
      keepRows.push_back(j);
    }
  }
  
  newStates = keepRows;
  NumericVector pi_w(p - NumZeros);
  NumericMatrix T_w(p - NumZeros, p - NumZeros);
  
  if (NumZeros == 0) {
    NumericMatrix diagonal(p,p);
    
    for (int i{0}; i < p; ++i) {
      diagonal(i,i) = 1.0 / Rw(i,0);
    }
    diagonal = matrix_product(diagonal, T);
    
    pi_w = pi;
    T_w = diagonal;
    
  }
  else {
    long n1{deleteRows.size()};
    long n2{keepRows.size()};
    
    NumericMatrix Spp(n2,n2);
    NumericMatrix Sp0(n2,n1);
    NumericMatrix S0p(n1,n2);
    NumericMatrix S00(n1,n1);
    
    NumericMatrix Taux(n2,n2);
    NumericMatrix diagonal(n2,n2);
    
    NumericMatrix pi0(1,n1);
    NumericMatrix pip(1,n2);
    
    NumericMatrix piaux(1,n2);
    
    for (int i{0}; i < n2; i++) {
      for (int j = 0; j < n2; j++) {
        Spp(i,j) = T(keepRows[i],keepRows[j]);
      }
      for (int j{0}; j < n1; j++) {
        Sp0(i,j) = T(keepRows[i],deleteRows[j]);
      }
      pip(0,i) = m_pi(0,keepRows[i]);
    }
    for (int i{0}; i < n1; i++) {
      for (int j{0}; j < n2; j++) {
        S0p(i,j) = T(deleteRows[i],keepRows[j]);
      }
      for (int j{0}; j < n1; j++){
        S00(i,j) = T(deleteRows[i],deleteRows[j]);
      }
      pi0(0,i) = m_pi(0,deleteRows[i]);
    }
    
    piaux = matrix_sum(pip, matrix_product(pi0, matrix_product(matrix_inverse(S00 * (-1.0)), S0p)));
    
    Taux = matrix_sum(Spp, matrix_product(Sp0, matrix_product(matrix_inverse(S00 * (-1.0)), S0p)));
    
    for (int i{0}; i < n2; ++i) {
      diagonal(i,i) = 1.0 / Rw(keepRows[i],0);
      pi_w[i] = piaux(0, i);
    }
    diagonal = matrix_product(diagonal, Taux);
    
    T_w = diagonal;
    
  }
  
  List L = List::create(Named("piw") = pi_w, _["Tw"] = T_w, _["new_states"] = keepRows);
  
  return L;
}


//' Second EM in the algorithm of bruer
// [[Rcpp::export]]
void secondEMstep(const NumericMatrix & observations, const NumericVector & weight, const NumericMatrix & censored, const NumericVector & rcweight, NumericVector & pi, NumericMatrix & T, NumericMatrix & R) {
  double lowerbound{1.0E-15}; //Makes a reward zero if it is below of this value - Otherwise it will never be zero
  long p{T.nrow()};
  long dim{R.ncol()};
  
  NumericMatrix Zmean(p,dim);
  NumericMatrix Zsum(p,1);
  
  NumericVector w(dim);
  
  //E-step
  for (int j{0}; j < dim; ++j) {
    
    w = w *0;
    w[j] = 1;
    
    List PHmarginal(linear_combination(w, pi, T, R));
    
    NumericVector m_piw = PHmarginal["piw"];
    
    NumericMatrix piw(1,m_piw.size(), m_piw.begin());
    
    NumericMatrix Tw = PHmarginal["Tw"];
    
    IntegerVector newStates = PHmarginal["new_states"];
    
    long pmarginal{Tw.nrow()};
    
    NumericVector m_e(pmarginal, 1);
    NumericMatrix e(pmarginal, 1, m_e.begin());
    
    NumericMatrix tw = matrix_product(Tw * (-1), e);
    
    NumericMatrix bvector(pmarginal, 1);
    NumericMatrix cmatrix(pmarginal, pmarginal);
    NumericMatrix aux_exp(pmarginal, pmarginal);
    
    NumericMatrix J(2 * pmarginal, 2 * pmarginal);
    NumericMatrix tProductPi(pmarginal, pmarginal);
    tProductPi = matrix_product(tw, piw);
    
    double density{0.0};
    
    //  Unccensored data
    for (int k{0}; k < observations.nrow(); ++k) {
      J = matrix_exponential(matrix_VanLoan(Tw, Tw, tProductPi) * observations(k, j));
      
      for (int i{0}; i < pmarginal; ++i) {
        for (int j{0}; j < pmarginal; ++j) {
          aux_exp(i,j) = J(i,j);
          cmatrix(i,j) = J(i,j + pmarginal);
        }
      }
      
      bvector = matrix_product(aux_exp, tw);
      density = matrix_product(piw, bvector)(0,0);
      
      //E-step
      for (int i{0}; i < pmarginal; ++i) {
        Zmean(newStates[i],j) += cmatrix(i,i) * weight[k] / density;
      }
    }
    
    
    //  Right-Censored Data
    if (censored.ncol() > 0) {
      tProductPi = matrix_product(e, piw);
    }
    for (int k{0}; k < censored.ncol(); ++k) {
      
      J = matrix_exponential(matrix_VanLoan(T, T, tProductPi) * censored(k,j));
      
      for (int i{0}; i < pmarginal; ++i) {
        for (int j{0}; j < pmarginal; ++j) {
          aux_exp(i,j) = J(i,j);
          cmatrix(i,j) = J(i,j + pmarginal);
        }
      }
      
      bvector = matrix_product(aux_exp, e);
      density = matrix_product(piw, bvector)(0, 0);
      
      //E-step
      for (int i{0}; i < pmarginal; ++i) {
        Zmean(newStates[i],j) += cmatrix(i,i) * rcweight[k] / density;
      }
    }
  }
  
  //M-step
  for (int i{0}; i < p; ++i) {
    for (int j{0}; j < dim; ++j) {
      Zsum(i,0) += Zmean(i,j);
    }
  }
  
  for (int j{0}; j < dim; ++j) {
    for (int i{0}; i < p; ++i){
      R(i,j) = Zmean(i,j) / Zsum(i,0);
      if (R(i,j) < lowerbound) {
        R(i,j) = 0;
      }
    }
  }
}

//' Sum data for input in Breur
// [[Rcpp::export]]
NumericVector sum_data(NumericMatrix x) {
  NumericVector sum_x(x.nrow());
  for (int i{0}; i < x.nrow(); ++i) {
    for (int j{0}; j < x.ncol(); ++j)
      sum_x[i] += x(i,j);
  }
  return sum_x;
}


////////////////////////////////////////////
// SOME SCALED VERSIONS (for regression):
////////////////////////////////////////////


//' Loglikelihood using RK, with scale
// [[Rcpp::export]]
double logLikelihoodPH_RKs(double h, NumericVector & pi, NumericMatrix & T, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * obs[0];
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) + log(scale1[k]));
    if (k < obs.size() - 1){dt = scale1[k + 1] * obs[k + 1] - scale1[k] * obs[k];}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * rcens[0];
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * rcens[k + 1] - scale2[k] * rcens[k];}
  }
  
  return logLh;
}



//' Loglikelihood using RK, with scale
// [[Rcpp::export]]
double logLikelihoodPH_RKs2(double h, NumericVector & pi, NumericMatrix & T, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * obs[0];
  }
  for (int k{0}; k < obs.size(); ++k) {
    density = matrix_product(m_pi, matrix_product(matrix_exponential(T * scale1[k] * obs[k]), t))(0,0);
    logLh += weight[k] * (log(density) + log(scale1[k]));
    //dt = scale1[k + 1] * (obs[k + 1] - obs[k]);
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * rcens[0];
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    //a_rungekutta(avector, dt, h, T);
    density = matrix_product(m_pi, matrix_product(matrix_exponential(T * scale1[k] * rcens[k]), e))(0,0);
    logLh += rcweight[k] * log(density);
    //dt = scale2[k + 1] * (rcens[k + 1] - rcens[k]);
  }
  
  return logLh;
}


//' Loglikelihood of matrix Weibull using RK
//' This is the fastest option
// [[Rcpp::export]]
double logLikelihoodMWeib_RKs(double h, NumericVector & pi, NumericMatrix & T, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * pow(obs[0], beta);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) + log(scale1[k]) + log(beta) + (beta -1) * log(obs[k]));
    if (k < obs.size() - 1){dt = scale1[k + 1] * pow(obs[k + 1], beta) - scale1[k] * pow(obs[k], beta);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * pow(rcens[0], beta);
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * pow(rcens[k + 1], beta) - scale2[k] * pow(rcens[k], beta);}
  }
  
  return logLh;
}



//' Loglikelihood of matrix Pareto using RK
// [[Rcpp::export]]
double logLikelihoodMPar_RKs(double h, NumericVector & pi, NumericMatrix & T, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * log(obs[0] / beta + 1);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) + log(scale1[k]) - log(obs[k] + beta));
    if (k < obs.size() - 1){dt = scale1[k + 1] * log(obs[k + 1] / beta + 1) - scale1[k] * log(obs[k] / beta + 1);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * log(rcens[0] / beta + 1);
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * log(rcens[k + 1] / beta + 1) - scale2[k] * log(rcens[k] / beta + 1);}
  }
  
  return logLh;
}


//' Loglikelihood of matrix Log-Logistic using RK
// [[Rcpp::export]]
double logLikelihoodMLogLogistic_RKs(double h, NumericVector & pi, NumericMatrix & T, NumericVector beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * log(pow(obs[0] / beta[0], beta[1]) + 1);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) + log(scale1[k]) + log(beta[1]) - log(beta[0]) + (beta[1] - 1) * (log(obs[k]) - log(beta[0])) - log(pow(obs[k] / beta[0], beta[1]) + 1));
    if (k < obs.size() - 1){dt = scale1[k + 1] * log(pow(obs[k + 1] / beta[0], beta[1]) + 1) - scale1[k] * log(pow(obs[k] / beta[0], beta[1]) + 1);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * log(pow(rcens[0] / beta[0], beta[1]) + 1);
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * log(pow(rcens[k + 1] / beta[0], beta[1]) + 1) - scale2[k] * log(pow(rcens[k] / beta[0], beta[1]) + 1);}
  }
  
  return logLh;
}

//' Loglikelihood of matrix Gompertz using RK
// [[Rcpp::export]]
double logLikelihoodMGomp_RKs(double h, NumericVector & pi, NumericMatrix & T, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * (exp(obs[0] * beta) - 1) / beta;
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) + log(scale1[k]) + obs[k] * beta);
    if (k < obs.size() - 1){dt = scale1[k + 1] * (exp(obs[k + 1] * beta) - 1) / beta - scale1[k] * (exp(obs[k] * beta) - 1) / beta;}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * (exp(rcens[0] * beta) - 1) / beta;
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * (exp(rcens[k + 1] * beta) - 1) / beta - scale2[k] * (exp(rcens[k] * beta) - 1) / beta;}
  }
  
  return logLh;
}


//' Loglikelihood of matrix GEV using RK
//' I am assuming that the sample is given in an increasing order
// [[Rcpp::export]]
double logLikelihoodMGEV_RKs(double h, NumericVector & pi, NumericMatrix & T, NumericVector beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  long N{obs.size()};
  
  if (beta[2] == 0) {
    // Non censored data
    if (N > 0) {
      dt = scale1[N - 1] * exp(-(obs[N - 1] - beta[0]) / beta[1]);
    }
    for (int k{1}; k <= N; ++k) {
      if(dt > 0) a_rungekutta(avector, dt, h, T);
      density = matrix_product(avector, t)(0,0);
      logLh += weight[N - k] * (log(density) + log(scale1[N - k]) - log(beta[1]) - (obs[N - k] - beta[0]) / beta[1]);
      if (k < N){dt = scale1[N - k - 1] * exp(-(obs[N - k - 1] - beta[0]) / beta[1]) - scale1[N - k] * exp(-(obs[N - k] - beta[0]) / beta[1]);}
    }
    //Right censored data
    N = rcens.size();
    if (N > 0) {
      dt = scale2[N - 1] * exp(-(rcens[N - 1] - beta[0]) / beta[1]);
      avector = clone(m_pi);
    }
    for (int k{1}; k <= N; ++k) {
      if(dt > 0) a_rungekutta(avector, dt, h, T);
      density = matrix_product(avector, e)(0,0);
      logLh += rcweight[N - k] * log(density);
      if (k < N){dt = scale2[N - k - 1] * exp(-(rcens[N - k - 1] - beta[0]) / beta[1]) - scale2[N - k] * exp(-(rcens[N - k] - beta[0]) / beta[1]);}
    }
  }
  else {
    // Non censored data
    if (N > 0) {
      dt = scale1[N - 1] * pow(1 + (beta[2] / beta[1]) * (obs[N - 1] - beta[0]) , - 1 / beta[2]);
    }
    for (int k{1}; k <= N; ++k) {
      if(dt > 0) a_rungekutta(avector, dt, h, T);
      density = matrix_product(avector, t)(0,0);
      logLh += weight[N - k] * (log(density) + log(scale1[N - k]) - log(beta[1]) - (1 + 1 / beta[2]) * log(1 + (beta[2] / beta[1]) * (obs[N - k] - beta[0])));
      if (k < N){dt = scale1[N - k - 1] * pow(1 + (beta[2] / beta[1]) * (obs[N - k - 1] - beta[0]) , - 1 / beta[2]) - scale1[N - k] * pow(1 + (beta[2] / beta[1]) * (obs[N - k] - beta[0]) , - 1 / beta[2]);}
    }
    //Right censored data
    N = rcens.size();
    if (N > 0) {
      dt = scale2[N - 1] * pow(1 + (beta[2] / beta[1]) * (rcens[N - 1] - beta[0]) , - 1 / beta[2]);
      avector = clone(m_pi);
    }
    for (int k{1}; k <= N; ++k) {
      if(dt > 0) a_rungekutta(avector, dt, h, T);
      density = matrix_product(avector, e)(0,0);
      logLh += rcweight[N - k] * log(density);
      if (k < N){dt = scale2[N - k - 1] * pow(1 + (beta[2] / beta[1]) * (rcens[N - k - 1] - beta[0]) , - 1 / beta[2]) - scale2[N - k] * pow(1 + (beta[2] / beta[1]) * (rcens[N - k] - beta[0]) , - 1 / beta[2]);}
    }
  }
  return logLh;
}

////////////////////////////////////////////
// SOME SCALED VERSIONS (for regression (AFT)):
////////////////////////////////////////////


//' Loglikelihood using RK, with scale
// [[Rcpp::export]]
double logLikelihoodPH_RKs1(double h, NumericVector & pi, NumericMatrix & T, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * obs[0];
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) + log(scale1[k]));
    if (k < obs.size() - 1){dt = scale1[k + 1] * obs[k + 1] - scale1[k] * obs[k];}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * rcens[0];
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * rcens[k + 1] - scale2[k] * rcens[k];}
  }
  
  return logLh;
}



//' Loglikelihood of matrix Weibull using RK
//' This is the fastest option
// [[Rcpp::export]]
double logLikelihoodMWeib_RKs1(double h, NumericVector & pi, NumericMatrix & T, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = pow(scale1[0] * obs[0], beta);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) + log(scale1[k]) + log(beta) + (beta -1) * log(scale1[k] * obs[k]));
    if (k < obs.size() - 1){dt = pow(scale1[k + 1] * obs[k + 1], beta) - pow(scale1[k] * obs[k], beta);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = pow(scale2[0] * rcens[0], beta);
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1){dt = pow(scale2[k + 1] * rcens[k + 1], beta) - pow(scale2[k] * rcens[k], beta);}
  }
  
  return logLh;
}



//' Loglikelihood of matrix Pareto using RK
// [[Rcpp::export]]
double logLikelihoodMPar_RKs1(double h, NumericVector & pi, NumericMatrix & T, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = log(scale1[0] * obs[0] / beta + 1);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) + log(scale1[k]) - log(scale1[k] * obs[k] + beta));
    if (k < obs.size() - 1){dt = log(scale1[k + 1] * obs[k + 1] / beta + 1) - log(scale1[k] * obs[k] / beta + 1);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = log(scale2[0] * rcens[0] / beta + 1);
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1){dt = log(scale2[k + 1] * rcens[k + 1] / beta + 1) - log(scale2[k] * rcens[k] / beta + 1);}
  }
  
  return logLh;
}


//' Loglikelihood of matrix Log-Logistic using RK
// [[Rcpp::export]]
double logLikelihoodMLogLogistic_RKs1(double h, NumericVector & pi, NumericMatrix & T, NumericVector beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = log(pow(scale1[0] * obs[0] / beta[0], beta[1]) + 1);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) + log(scale1[k]) + log(beta[1]) - log(beta[0]) + (beta[1] - 1) * (log(scale1[k] * obs[k]) - log(beta[0])) - log(pow(scale1[k] * obs[k] / beta[0], beta[1]) + 1));
    if (k < obs.size() - 1){dt = log(pow(scale1[k + 1] * obs[k + 1] / beta[0], beta[1]) + 1) - log(pow(scale1[k] * obs[k] / beta[0], beta[1]) + 1);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = log(pow(scale2[0] * rcens[0] / beta[0], beta[1]) + 1);
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1){dt = log(pow(scale2[k + 1] * rcens[k + 1] / beta[0], beta[1]) + 1) - log(pow(scale2[k] * rcens[k] / beta[0], beta[1]) + 1);}
  }
  
  return logLh;
}

//' Loglikelihood of matrix Gompertz using RK
// [[Rcpp::export]]
double logLikelihoodMGomp_RKs1(double h, NumericVector & pi, NumericMatrix & T, double beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = (exp(scale1[0] * obs[0] * beta) - 1) / beta;
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) + log(scale1[k]) + scale1[k] * obs[k] * beta);
    if (k < obs.size() - 1){dt = (exp(scale1[k + 1] * obs[k + 1] * beta) - 1) / beta - (exp(scale1[k] * obs[k] * beta) - 1) / beta;}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = (exp(scale2[0] * rcens[0] * beta) - 1) / beta;
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1){dt = (exp(scale2[k + 1] * rcens[k + 1] * beta) - 1) / beta - (exp(scale2[k] * rcens[k] * beta) - 1) / beta;}
  }
  
  return logLh;
}


//' Loglikelihood of matrix GEV using RK
//' I am assuming that the sample is given in an increasing order
// [[Rcpp::export]]
double logLikelihoodMGEV_RKs1(double h, NumericVector & pi, NumericMatrix & T, NumericVector beta, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  long N{obs.size()};
  
  if (beta[2] == 0) {
    // Non censored data
    if (N > 0) {
      dt = exp(-(scale1[N - 1] * obs[N - 1] - beta[0]) / beta[1]);
    }
    for (int k{1}; k <= N; ++k) {
      if(dt > 0) a_rungekutta(avector, dt, h, T);
      density = matrix_product(avector, t)(0,0);
      logLh += weight[N - k] * (log(density) + log(scale1[N - k]) - log(beta[1]) - (scale1[N - k] * obs[N - k] - beta[0]) / beta[1]);
      if (k < N){dt = exp(-(scale1[N - k - 1] * obs[N - k - 1] - beta[0]) / beta[1]) - exp(-(scale1[N - k] * obs[N - k] - beta[0]) / beta[1]);}
    }
    //Right censored data
    N = rcens.size();
    if (N > 0) {
      dt = exp(-(scale2[N - 1] * rcens[N - 1] - beta[0]) / beta[1]);
      avector = clone(m_pi);
    }
    for (int k{1}; k <= N; ++k) {
      if(dt > 0) a_rungekutta(avector, dt, h, T);
      density = matrix_product(avector, e)(0,0);
      logLh += rcweight[N - k] * log(density);
      if (k < N){dt = exp(-(scale2[N - k - 1] * rcens[N - k - 1] - beta[0]) / beta[1]) - exp(-(scale2[N - k] * rcens[N - k] - beta[0]) / beta[1]);}
    }
  }
  else {
    // Non censored data
    if (N > 0) {
      dt = pow(1 + (beta[2] / beta[1]) * (scale1[N - 1] * obs[N - 1] - beta[0]) , - 1 / beta[2]);
    }
    for (int k{1}; k <= N; ++k) {
      if(dt > 0) a_rungekutta(avector, dt, h, T);
      density = matrix_product(avector, t)(0,0);
      logLh += weight[N - k] * (log(density) + log(scale1[N - k]) - log(beta[1]) - (1 + 1 / beta[2]) * log(1 + (beta[2] / beta[1]) * (scale1[N - k] * obs[N - k] - beta[0])));
      if (k < N){dt = pow(1 + (beta[2] / beta[1]) * (scale1[N - k - 1] * obs[N - k - 1] - beta[0]) , - 1 / beta[2]) - pow(1 + (beta[2] / beta[1]) * (scale1[N - k] * obs[N - k] - beta[0]) , - 1 / beta[2]);}
    }
    //Right censored data
    N = rcens.size();
    if (N > 0) {
      dt = pow(1 + (beta[2] / beta[1]) * (scale2[N - 1] * rcens[N - 1] - beta[0]) , - 1 / beta[2]);
      avector = clone(m_pi);
    }
    for (int k{1}; k <= N; ++k) {
      if(dt > 0) a_rungekutta(avector, dt, h, T);
      density = matrix_product(avector, e)(0,0);
      logLh += rcweight[N - k] * log(density);
      if (k < N){dt = pow(1 + (beta[2] / beta[1]) * (scale2[N - k - 1] * rcens[N - k - 1] - beta[0]) , - 1 / beta[2]) - pow(1 + (beta[2] / beta[1]) * (scale2[N - k] * rcens[N - k] - beta[0]) , - 1 / beta[2]);}
    }
  }
  return logLh;
}


////////////////////////////////////////////
// SOME SCALED VERSIONS (for double regression):
////////////////////////////////////////////


//' Loglikelihood of matrix Weibull using RK
//' This is the fastest option
// [[Rcpp::export]]
double logLikelihoodMWeib_RKs_double(double h, NumericVector & pi, NumericMatrix & T, const NumericVector & beta1, const NumericVector & beta2, const NumericVector & obs, const NumericVector & weight, const NumericVector & rcens, const NumericVector & rcweight, const NumericVector & scale1, const NumericVector & scale2) {
  long p{T.nrow()};
  NumericMatrix m_pi(1,p, pi.begin());
  
  NumericMatrix avector(1,p);
  
  NumericVector m_e(p, 1);
  NumericMatrix e(p, 1, m_e.begin());
  
  NumericMatrix t = matrix_product(T * (-1), e);
  
  // Uncensored data
  //   initial condition
  avector = clone(m_pi);
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * pow(obs[0], beta1[0]);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, t)(0,0);
    logLh += weight[k] * (log(density) + log(scale1[k]) + log(beta1[k]) + (beta1[k] -1) * log(obs[k]));
    if (k < obs.size() - 1){dt = scale1[k + 1] * pow(obs[k + 1], beta1[k + 1]) - scale1[k] * pow(obs[k], beta1[k]);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * pow(rcens[0], beta2[0]);
    avector = clone(m_pi);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, T);
    density = matrix_product(avector, e)(0,0);
    logLh += rcweight[k] * log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * pow(rcens[k + 1], beta2[k + 1]) - scale2[k] * pow(rcens[k], beta2[k]);}
  }
  
  return logLh;
}
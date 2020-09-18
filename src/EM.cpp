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
    
    dt = obs[k + 1] - obs[k];
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
    
    dt = rcens[k + 1] - rcens[k];
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
    dt = obs[k + 1] - obs[k];
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
    dt = rcens[k + 1] - rcens[k];
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
    g_inv_val2 = g_inv(obs[k + 1], beta);
    dt = g_inv_val2[0] - g_inv_val[0];
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
    g_inv_val2 = g_inv(rcens[k + 1], beta);
    dt = g_inv_val2[0] - g_inv_val[0];
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
    dt = pow(obs[k + 1], beta) - pow(obs[k], beta);
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
    dt = pow(rcens[k + 1], beta) - pow(rcens[k], beta);
  }
  
  return logLh;
}



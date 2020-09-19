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
    dt = log(obs[k + 1] / beta + 1) - log(obs[k] / beta + 1);
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
    dt = log(rcens[k + 1] / beta + 1) - log(rcens[k] / beta + 1);
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
    dt = (exp(obs[k + 1] * beta) - 1) / beta - (exp(obs[k] * beta) - 1) / beta;
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
    dt = (exp(rcens[k + 1] * beta) - 1) / beta - (exp(rcens[k] * beta) - 1) / beta;
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
      dt = exp(-(obs[N - k - 1] - beta[0]) / beta[1]) - exp(-(obs[N - k] - beta[0]) / beta[1]);
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
      dt = exp(-(rcens[N - k - 1] - beta[0]) / beta[1]) - exp(-(rcens[N - k] - beta[0]) / beta[1]);
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
      dt = pow(1 + (beta[2] / beta[1]) * (obs[N - k - 1] - beta[0]) , - 1 / beta[2]) - pow(1 + (beta[2] / beta[1]) * (obs[N - k] - beta[0]) , - 1 / beta[2]);
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
      dt = pow(1 + (beta[2] / beta[1]) * (rcens[N - k - 1] - beta[0]) , - 1 / beta[2]) - pow(1 + (beta[2] / beta[1]) * (rcens[N - k] - beta[0]) , - 1 / beta[2]);
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
  
  NumericMatrix m_alpha(1,p, alpha.begin());
  
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



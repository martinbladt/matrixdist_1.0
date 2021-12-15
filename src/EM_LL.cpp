# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]


//' Runge Kutta for the calculation of the a,b and c vectors in a EM step
//' 
//' Performs the RK of forth order.
//' 
//' @param avector The a vector.
//' @param bvector The b vector.
//' @param cmatrix The c matrix.
//' @param dt The increment.
//' @param h Step-length.
//' @param S Sub-intensity.
//' @param s Exit rates.
//' 
// [[Rcpp::export]]
void runge_kutta(arma::vec & avector, arma::mat & bvector, arma::mat & cmatrix, double dt, double h, const arma::mat & S, const arma::mat & s) {
  unsigned p{S.n_rows};
  int j{};
  int m{};
  double eps{};
  double sum{};
  
  int i{};
  i = dt / h;
  
  double h2{};
  h2 = dt / (i + 1);
  
  arma::mat ka(4,p);
  arma::mat kb(4,p);
  arma::mat kc1(p,p);
  arma::mat kc2(p,p);
  arma::mat kc3(p,p);
  arma::mat kc4(p,p);
  
  for (eps = 0; eps <= dt - h2 / 2; eps += h2) {
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(j,i) * avector[j];
      }
      ka(0,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(j,i) * (avector[j] + ka(0,j) / 2);
      }
      ka(1,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(j,i) * (avector[j] + ka(1,j) / 2);
      }
      ka(2,i) = h2 * sum;
    }
    for (i=0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j)
      {
        sum += S(j,i) * (avector[j] + ka(2,j));
      }
      ka(3,i) = h2 * sum;
    }
    
    
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(i,j) * bvector(j,0) ;
      }
      kb(0,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(i,j) * (bvector(j,0)  + kb(0,j) / 2);
      }
      kb(1,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(i,j) * (bvector(j,0)  + kb(1,j) / 2);
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
        sum = s(m,0) * avector[i];
        for (j = 0; j < p; ++j) {
          sum += S(m,j) * cmatrix(j,i);
        }
        kc1(m,i) = h2 * sum;
      }
    }
    for (m = 0; m < p; ++m) {
      for (i = 0; i < p; ++i) {
        sum = s(m,0) * (avector[i] + ka(0,i) / 2);
        for (j = 0; j < p; ++j) {
          sum += S(m,j) * (cmatrix(j,i) + kc1(j,i) / 2);
        }
        kc2(m,i) = h2 * sum;
      }
    }
    for (m = 0; m < p; ++m) {
      for (i = 0; i < p; ++i) {
        sum = s(m,0) * (avector[i] + ka(1,i) / 2);
        for (j = 0; j < p; ++j) {
          sum += S(m,j) * (cmatrix(j,i) + kc2(j,i) / 2);
        }
        kc3(m,i) = h2 * sum;
      }
    }
    for (m = 0; m < p; ++m) {
      for (i = 0; i < p; ++i) {
        sum = s(m,0) * (avector[i] + ka(2,i));
        for (j = 0; j < p; ++j) {
          sum += S(m,j) * (cmatrix(j,i) + kc3(j,i));
        }
        kc4(m,i) = h2 * sum;
      }
    }
    
    for (i = 0; i < p; ++i) {
      avector[i] += (ka(0,i) + 2 * ka(1,i) + 2 * ka(2,i) + ka(3,i)) / 6;
      bvector(i,0) += (kb(0,i) + 2 * kb(1,i) + 2 * kb(2,i) + kb(3,i)) / 6;
      for (j = 0; j < p; ++j) {
        cmatrix(i,j) +=(kc1(i,j) + 2 * kc2(i,j) + 2 * kc3(i,j) + kc4(i,j)) / 6;
      }
    }
  }
}


//' EM step using Runge Kutta
//' 
//' Computes one step of the EM algorithm by using a Runge-Kutta method of 4th order.
//' 
//' @param h Step-length.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param obs The observations.
//' @param weight The weights for the observations.
//' @param rcens Censored observations.
//' @param rcweight The weights for the censored observations.
//' 
// [[Rcpp::export]]
void EMstep_RK(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  unsigned p{S.n_rows};
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat Bmean = arma::zeros(p,1);
  arma::mat Zmean = arma::zeros(p,1);
  arma::mat Nmean = arma::zeros(p,p + 1);
  
  arma::vec avector(p);
  arma::mat bvector(p,1);
  arma::mat cmatrix(p,p);
  arma::mat aux_exp(p,p);
  
  arma::mat aux_mat(1,1);
  
  // initial conditions
  avector = alpha;
  bvector = t;
  
  double dt{0.0};
  if (obs.size() > 0) {
    dt = obs[0];
  }
  
  double SumOfWeights{0.0};
  double density{0.0};
  
  // E step
  //  Uncensored data
  for (int k{0}; k < obs.size(); ++k) {
    SumOfWeights += weight[k];
    
    runge_kutta(avector, bvector, cmatrix, dt, h, S, t);
    aux_mat = alpha.t() * bvector;
    density = aux_mat(0,0);
    
    // E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha[i] * bvector(i,0) * weight[k] / density;
      Nmean(i,p) += avector[i] * t(i,0) * weight[k] / density;
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
    cmatrix = cmatrix * 0;
    avector = alpha;
    bvector = e;
  }
  for (int k{0}; k < rcens.size(); ++k) {
    SumOfCensored += rcweight[k];
    
    runge_kutta(avector, bvector, cmatrix, dt, h, S, e);
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
//' @param avector The a vector.
//' @param dt Increment.
//' @param h Step-length.
//' @param S Sub-intensity.
//' 
// [[Rcpp::export]]
void a_rungekutta(arma::vec & avector, double dt, double h, const arma::mat & S) {
  unsigned p{S.n_rows};
  int j{};
  double eps{};
  double sum{};
  
  int i{};
  i = dt / h;
  
  double h2{};
  h2 = dt / (i + 1);
  
  arma::mat ka(4,p);
  
  for (eps = 0; eps <= dt - h2 / 2; eps += h2) {
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(j,i) * avector[j];
      }
      ka(0,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(j,i) * (avector[j] + ka(0,j) / 2);
      }
      ka(1,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(j,i) * (avector[j] + ka(1,j) / 2);
      }
      ka(2,i) = h2 * sum;
    }
    for (i = 0; i < p; ++i) {
      sum = 0;
      for (j = 0; j < p; ++j) {
        sum += S(j,i) * (avector[j] + ka(2,j));
      }
      ka(3,i) = h2 * sum;
    }
    
    for (i = 0; i < p; ++i) {
      avector[i] += (ka(0,i) + 2 * ka(1,i) + 2 * ka(2,i) + ka(3,i)) / 6;
    }
  }
}


//' Loglikelihood using RK
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Step-length.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param obs The observations.
//' @param weight Weight of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weight of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodPH_RK(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  unsigned p{S.n_rows};
 
  arma::vec avector(p);
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  // Initial condition
  avector = alpha;
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = obs[0];
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * t;
    density = aux_mat(0,0);
    logLh += weight[k] * std::log(density);
    if (k < obs.size() - 1) {dt = obs[k + 1] - obs[k]; }
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = rcens[0];
    avector = alpha;
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1) {dt = rcens[k + 1] - rcens[k];}
  }
  
  return logLh;
}


//' Loglikelihood of matrix Weibull using RK
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Step-length.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weight of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weight of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMweibull_RK(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  
  arma::vec avector(p);
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  // Initial condition
  avector = alpha;
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = pow(obs[0], beta);
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta) + (beta -1) * std::log(obs[k]));
    if (k < obs.size() - 1) {dt = pow(obs[k + 1], beta) - pow(obs[k], beta);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = pow(rcens[0], beta);
    avector = alpha;
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1) {dt = pow(rcens[k + 1], beta) - pow(rcens[k], beta);}
  }
  
  return logLh;
}


//' Loglikelihood of matrix Pareto using RK
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Step-length.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weight of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weight of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMpareto_RK(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  
  arma::vec avector(p);
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);

  // Initial condition
  avector = alpha;
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = std::log(obs[0] / beta + 1);
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) - std::log(obs[k] + beta));
    if (k < obs.size() - 1) {dt = std::log(obs[k + 1] / beta + 1) - std::log(obs[k] / beta + 1);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = std::log(rcens[0] / beta + 1);
    avector = alpha;
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1) {dt = std::log(rcens[k + 1] / beta + 1) - std::log(rcens[k] / beta + 1);}
  }
  
  return logLh;
}


//' Loglikelihood of matrix Lognormal using RK
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Step-length.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weight of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weight of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMlognormal_RK(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  
  arma::vec avector(p);
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  // Initial condition
  avector = alpha;
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = pow(std::log(obs[0] + 1), beta);
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta) + (beta -1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
    if (k < obs.size() - 1) {dt = pow(std::log(obs[k + 1] + 1), beta) - pow(std::log(obs[k] + 1), beta);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = pow(std::log(rcens[0] + 1), beta);
    avector = alpha;
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1) {dt = pow(std::log(rcens[k + 1] + 1), beta) - pow(std::log(rcens[k] + 1), beta);}
  }
  
  return logLh;
}


//' Loglikelihood of matrix Log-logistic using RK
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Step-length.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param beta Parameters of transformation.
//' @param obs The observations.
//' @param weight Weight of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weight of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMloglogistic_RK(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta[0] < 0 || beta[1] < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  
  arma::vec avector(p);
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  // Initial condition
  avector = alpha;
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = std::log(pow(obs[0] / beta[0], beta[1]) + 1);
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta[1]) - std::log(beta[0]) + (beta[1] - 1) * (std::log(obs[k]) - std::log(beta[0])) - std::log(pow(obs[k] / beta[0], beta[1]) + 1));
    if (k < obs.size() - 1) {dt = std::log(pow(obs[k + 1] / beta[0], beta[1]) + 1) - std::log(pow(obs[k] / beta[0], beta[1]) + 1);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = std::log(pow(rcens[0] / beta[0], beta[1]) + 1);
    avector = alpha;
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1) {dt = std::log(pow(rcens[k + 1] / beta[0], beta[1]) + 1) - std::log(pow(rcens[k] / beta[0], beta[1]) + 1);}
  }
  
  return logLh;
}


//' Loglikelihood of matrix Gompertz using RK
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Step-length.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weight of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weight of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMgompertz_RK(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  
  arma::vec avector(p);
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  // Initial condition
  avector = alpha;
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = (exp(obs[0] * beta) - 1) / beta;
  }
  for (int k{0}; k < obs.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + obs[k] * beta);
    if (k < obs.size() - 1) {dt = (exp(obs[k + 1] * beta) - 1) / beta - (exp(obs[k] * beta) - 1) / beta;}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = (exp(rcens[0] * beta) - 1) / beta;
    avector = alpha;
  }
  for (int k{0}; k < rcens.size(); ++k) {
    a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1) {dt = (exp(rcens[k + 1] * beta) - 1) / beta - (exp(rcens[k] * beta) - 1) / beta;}
  }
  
  return logLh;
}


//' Loglikelihood of matrix GEV using RK
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Step-length.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param beta Parameter of transformation
//' @param obs The observations.
//' @param weight Weight of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weight of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMgev_RK(double h, arma::vec  & alpha, arma::mat & S, Rcpp::NumericVector beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta[1] < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  
  arma::vec avector(p);
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  // Initial condition
  avector = alpha;
  
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
      aux_mat = avector.t() * t;
      density = aux_mat(0,0);
      logLh += weight[N - k] * (std::log(density) - std::log(beta[1]) - (obs[N - k] - beta[0]) / beta[1]);
      if (k < N) {dt = exp(-(obs[N - k - 1] - beta[0]) / beta[1]) - exp(-(obs[N - k] - beta[0]) / beta[1]);}
    }
    //Right censored data
    N = rcens.size();
    if (N > 0) {
      dt = exp(-(rcens[N - 1] - beta[0]) / beta[1]);
      avector = alpha;
    }
    for (int k{1}; k <= N; ++k) {
      a_rungekutta(avector, dt, h, S);
      aux_mat = avector.t() * e;
      density = aux_mat(0,0);
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
      aux_mat = avector.t() * t;
      density = aux_mat(0,0);
      logLh += weight[N - k] * (std::log(density) - std::log(beta[1]) - (1 + 1 / beta[2]) * std::log(1 + (beta[2] / beta[1]) * (obs[N - k] - beta[0])));
      if (k < N) {dt = pow(1 + (beta[2] / beta[1]) * (obs[N - k - 1] - beta[0]) , - 1 / beta[2]) - pow(1 + (beta[2] / beta[1]) * (obs[N - k] - beta[0]) , - 1 / beta[2]);}
      
    }
    //Right censored data
    N = rcens.size();
    if (N > 0) {
      dt = pow(1 + (beta[2] / beta[1]) * (rcens[N - 1] - beta[0]) , - 1 / beta[2]);
      avector = alpha;
    }
    for (int k{1}; k <= N; ++k) {
      a_rungekutta(avector, dt, h, S);
      aux_mat = avector.t() * e;
      density = aux_mat(0,0);
      logLh += rcweight[N - k] * std::log(density);
      if (k < N) {dt = pow(1 + (beta[2] / beta[1]) * (rcens[N - k - 1] - beta[0]) , - 1 / beta[2]) - pow(1 + (beta[2] / beta[1]) * (rcens[N - k] - beta[0]) , - 1 / beta[2]);}
    }
  }
  return logLh;
}


////////////////////////////////////////////////////////
// Scaled versions of loglikelihoods (for regression)://
///////////////////////////////////////////////////////

//' Loglikelihood of PH using RK
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Step-length.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param obs The observations.
//' @param weight Weight of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weight of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodPH_RKs(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  unsigned p{S.n_rows};
  
  arma::vec avector(p);
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  // Initial condition
  avector = alpha;
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * obs[0];
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]));
    if (k < obs.size() - 1){dt = scale1[k + 1] * obs[k + 1] - scale1[k] * obs[k];}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * rcens[0];
    avector = alpha;
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * rcens[k + 1] - scale2[k] * rcens[k];}
  }
  
  return logLh;
}


//' Loglikelihood of matrix-Weibull using RK
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Step-length.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weight of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weight of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMweibull_RKs(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  unsigned p{S.n_rows};
  
  arma::vec avector(p);
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  // Initial condition
  avector = alpha;
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * pow(obs[0], beta);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta) + (beta -1) * std::log(obs[k]));
    if (k < obs.size() - 1){dt = scale1[k + 1] * pow(obs[k + 1], beta) - scale1[k] * pow(obs[k], beta);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * pow(rcens[0], beta);
    avector = alpha;
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * pow(rcens[k + 1], beta) - scale2[k] * pow(rcens[k], beta);}
  }
  
  return logLh;
}


//' Loglikelihood of matrix-Pareto using RK
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Step-length.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weight of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weight of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMpareto_RKs(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  unsigned p{S.n_rows};
  
  arma::vec avector(p);
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  // Initial condition
  avector = alpha;
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * std::log(obs[0] / beta + 1);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) - std::log(obs[k] + beta));
    if (k < obs.size() - 1){dt = scale1[k + 1] * std::log(obs[k + 1] / beta + 1) - scale1[k] * std::log(obs[k] / beta + 1);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * std::log(rcens[0] / beta + 1);
    avector = alpha;
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * std::log(rcens[k + 1] / beta + 1) - scale2[k] * std::log(rcens[k] / beta + 1);}
  }
  
  return logLh;
}


//' Loglikelihood of matrix Lognormal using RK
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Step-length.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weight of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weight of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMlognormal_RKs(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  unsigned p{S.n_rows};
  
  arma::vec avector(p);
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  // Initial condition
  avector = alpha;
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * pow(std::log(obs[0] + 1), beta);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta) + (beta -1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
    if (k < obs.size() - 1){dt = scale1[k + 1] * pow(std::log(obs[k + 1] + 1), beta) - scale1[k] * pow(std::log(obs[k] + 1), beta);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * pow(std::log(rcens[0] + 1), beta);
    avector = alpha;
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * pow(std::log(rcens[k + 1] + 1), beta) - scale2[k] * pow(std::log(rcens[k] + 1), beta);}
  }
  
  return logLh;
}


//' Loglikelihood of matrix-loglogistic using RK
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Step-length.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param beta Parameters of transformation.
//' @param obs The observations.
//' @param weight Weight of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weight of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMloglogistic_RKs(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  unsigned p{S.n_rows};
  
  arma::vec avector(p);
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  // Initial condition
  avector = alpha;
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * std::log(pow(obs[0] / beta[0], beta[1]) + 1);
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta[1]) - std::log(beta[0]) + (beta[1] - 1) * (std::log(obs[k]) - std::log(beta[0])) - std::log(pow(obs[k] / beta[0], beta[1]) + 1));
    if (k < obs.size() - 1){dt = scale1[k + 1] * std::log(pow(obs[k + 1] / beta[0], beta[1]) + 1) - scale1[k] * std::log(pow(obs[k] / beta[0], beta[1]) + 1);}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * std::log(pow(rcens[0] / beta[0], beta[1]) + 1);
    avector = alpha;
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * std::log(pow(rcens[k + 1] / beta[0], beta[1]) + 1) - scale2[k] * std::log(pow(rcens[k] / beta[0], beta[1]) + 1);}
  }
  
  return logLh;
}


//' Loglikelihood of matrix Gompertz using RK
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Step-length.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weight of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weight of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMgompertz_RKs(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  unsigned p{S.n_rows};
  
  arma::vec avector(p);
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  // Initial condition
  avector = alpha;
  
  double dt{0.0};
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  if (obs.size() > 0) {
    dt = scale1[0] * (exp(obs[0] * beta) - 1) / beta;
  }
  for (int k{0}; k < obs.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * t;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + obs[k] * beta);
    if (k < obs.size() - 1){dt = scale1[k + 1] * (exp(obs[k + 1] * beta) - 1) / beta - scale1[k] * (exp(obs[k] * beta) - 1) / beta;}
  }
  //Right censored data
  if (rcens.size() > 0) {
    dt = scale2[0] * (exp(rcens[0] * beta) - 1) / beta;
    avector = alpha;
  }
  for (int k{0}; k < rcens.size(); ++k) {
    if(dt > 0) a_rungekutta(avector, dt, h, S);
    aux_mat = avector.t() * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
    if (k < rcens.size() - 1){dt = scale2[k + 1] * (exp(rcens[k + 1] * beta) - 1) / beta - scale2[k] * (exp(rcens[k] * beta) - 1) / beta;}
  }
  
  return logLh;
}

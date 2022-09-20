#include "m_exp.h"
# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]

//' EM for discrete phase-type
//' 
//' @param alpha Initial probabilities.
//' @param S Sub-transition matrix.
//' @param obs The observations.
//' @param weight The weights for the observations.
//' 
// [[Rcpp::export]]
void EMstep_dph(arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight) {
  unsigned p{S.n_rows};
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = e - (S * e) ;
  
  arma::mat Bmean = arma::zeros(p,1);
  arma::mat Nmean = arma::zeros(p,p + 1);

  arma::mat avector(1,p);
  arma::mat bvector(p,1);
  arma::mat cmatrix(p,p);
  arma::mat aux_pow(p,p);
  
  arma::mat aux_mat(1,1);
  
  arma::mat J(2 * p,2 * p);
  arma::mat s_prod_alpha(p,p);
  s_prod_alpha = exit_vect * alpha.t();
  
  J = matrix_vanloan(S, S, s_prod_alpha);
  
  double max_val{max(obs)};
  std::vector<arma::mat> vect = vector_of_powers(J, max_val);
  
  double sum_weights{0.0};
  double density{0.0};
  
  //E-step
  for (int k{0}; k < obs.size(); ++k) {
    sum_weights += weight[k];
    
    //Matrix power
    J = vect[obs[k] - 1];
    
    //Separate matrix
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_pow(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }
    
    avector = alpha.t() * aux_pow;
    bvector = aux_pow * exit_vect;
    aux_mat = alpha.t() * bvector;
    density = aux_mat(0,0);
    
    //E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha[i] * bvector(i,0) * weight[k] / density;
      Nmean(i,p) += avector(0,i) * exit_vect(i,0) * weight[k] / density;
      if (obs[k] > 1) {
        for (int j{0}; j < p; ++j) {
          Nmean(i,j) += S(i,j) * cmatrix(j,i) * weight[k] / density;
        }
      }
    }
  }
  
  arma::colvec Ncum = arma::sum(Nmean, 1);
  
  // M step
  for (int i{0}; i < p; ++i) {
    alpha[i] = Bmean(i,0) / (sum_weights);
    for (int j{0}; j < p; ++j) {
      S(i,j) = Nmean(i,j) / Ncum[i];
    }
  }
}


//' EM for discrete phase-type MoE
//' 
//' @param alpha Initial probabilities.
//' @param S Sub-transition matrix.
//' @param obs The observations.
//' @param weight The weights for the observations.
//' 
// [[Rcpp::export]]
arma::mat EMstep_dph_MoE(arma::mat & alpha,  arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight) {
  unsigned p{S.n_rows};
  
  arma::mat Bmatrix(obs.size(), p);
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = e - (S * e) ;
  
  arma::rowvec alpha_aux(alpha.row(0));
  
  arma::mat Nmean = arma::zeros(p,p + 1);
  
  arma::mat avector(1,p);
  arma::mat bvector(p,1);
  arma::mat cmatrix(p,p);
  arma::mat aux_pow(p,p);
  
  arma::mat aux_mat(1,1);
  
  arma::mat s_prod_alpha(p,p);
  
  double max_val{max(obs)};
  std::vector<arma::mat> vect = vector_of_powers(S, max_val);
  
  double density{0.0};
  
  //E-step
  for (int k{0}; k < obs.size(); ++k) {
    alpha_aux = alpha.row(k);
    s_prod_alpha = exit_vect * alpha_aux;
    
    //Matrix power
    aux_pow = vect[obs[k] - 1];
    
    avector = alpha_aux * aux_pow;
    bvector = aux_pow * exit_vect;
    aux_mat = alpha_aux * bvector;
    density = aux_mat(0,0);
    
    //E-step
    for (int i{0}; i < p; ++i) {
      Bmatrix(k,i) += alpha_aux[i] * bvector(i,0) * weight[k] / density;
      Nmean(i,p) += avector(0,i) * exit_vect(i,0) * weight[k] / density;
      if (obs[k] > 1) {
        cmatrix = arma::zeros(p,p);
        for (int l{0}; l <= obs[k] - 2; ++l) {
          cmatrix = cmatrix + vect[obs[k] - 2 - l] * s_prod_alpha * vect[l];
        }
        for (int j{0}; j < p; ++j) {
          Nmean(i,j) += S(i,j) * cmatrix(j,i) * weight[k] / density;
        }
      }
    }
  }
  
  arma::colvec Ncum = arma::sum(Nmean, 1);
  
  // M step
  for (int i{0}; i < p; ++i) {
    for (int j{0}; j < p; ++j) {
      S(i,j) = Nmean(i,j) / Ncum[i];
    }
  }
  return(Bmatrix);
}



//' EM for discrete bivariate phase-type
//' 
//' @param alpha Initial probabilities.
//' @param S11 Sub-transition matrix.
//' @param S12 Matrix.
//' @param S22 Sub-transition matrix.
//' @param obs The observations.
//' @param weight The weights for the observations.
//' 
// [[Rcpp::export]]
void EMstep_bivdph(arma::vec & alpha, arma::mat & S11, arma::mat & S12, arma::mat & S22, const Rcpp::NumericMatrix & obs, const Rcpp::NumericVector & weight) {
  unsigned p1{S11.n_rows};
  unsigned p2{S22.n_rows};
  unsigned p{p1 + p2};
  
  arma::mat e;
  e.ones(S22.n_cols, 1);
  arma::mat exit_vect = e - (S22 * e);
  
  arma::mat Bmean = arma::zeros(p1,1);
  arma::mat Nmean = arma::zeros(p,p + 1);
  
  arma::mat aux_matrix0(p1,p2);
  arma::mat aux_matrix1(p1,1);
  arma::mat aux_matrix2(p2,p1);
  arma::mat aux_matrix3(1,p2);
  
  arma::mat cmatrix1(p1,p1);
  arma::mat cmatrix2(p2,p2);
  
  arma::mat aux_mat(1,1);
  
  double max_val1{max(obs.column(0))};
  double max_val2{max(obs.column(1))};
  
  std::vector<arma::mat> vect1 = vector_of_powers(S11, max_val1);
  std::vector<arma::mat> vect2 = vector_of_powers(S22, max_val2);
  
  double sum_weights{0.0};
  double density{0.0};
  
  //E-step
  for (int k{0}; k < obs.nrow(); ++k) {
    sum_weights += weight[k];
    
    aux_matrix0 = vect1[obs(k, 0) - 1] * S12 * vect2[obs(k, 1) - 1];
    aux_matrix1 = aux_matrix0 * exit_vect;
    aux_matrix2 = vect2[obs(k, 1) - 1] * exit_vect * alpha.t() * vect1[obs(k, 0) - 1];
    aux_matrix3 = alpha.t() * aux_matrix0;

    aux_mat = alpha.t() * aux_matrix1;
    density = aux_mat(0,0);
    
    //E-step
    for (int i{0}; i < p1; ++i) {
      Bmean(i,0) += alpha[i] * aux_matrix1(i,0) * weight[k] / density;
      if (obs(k, 0) > 1) {
        cmatrix1 = cmatrix1 * 0;
        for (int m{0}; m <= obs(k, 0) - 2; ++m) {
          cmatrix1 = cmatrix1 + vect1[obs(k, 0) - m - 2] * S12 * vect2[obs(k, 1) - 1] * exit_vect * alpha.t() * vect1[m];
        }
        for (int j{0}; j < p1; ++j) {
          Nmean(i,j) += S11(i,j) * cmatrix1(j,i) * weight[k] / density;
        }
      }
      for (int j{0}; j < p2; ++j) {
        Nmean(i,j + p1) += S12(i,j) * aux_matrix2(j, i) * weight[k] / density;
      }
    }
    
    for (int i{0}; i < p2; ++i) {
      Nmean(i + p1,p) += aux_matrix3(0,i) * exit_vect(i,0) * weight[k] / density;
      if (obs(k, 1) > 1) {
        cmatrix2 = cmatrix2 * 0;
        for (int m{0}; m <= obs(k, 1) - 2; ++m) {
          cmatrix2 = cmatrix2 + vect2[obs(k, 1) - m - 2] * exit_vect * alpha.t() * vect1[obs(k, 0) - 1] * S12 * vect2[m];
        }
        for (int j{0}; j < p2; ++j) {
          Nmean(i + p1,j + p1) += S22(i,j) * cmatrix2(j,i) * weight[k] / density;
        }
      }
    }
  }
  
  arma::colvec Ncum = arma::sum(Nmean, 1);
  
  // M step
  for (int i{0}; i < p1; ++i) {
    alpha[i] = Bmean(i,0) / (sum_weights);
    for (int j{0}; j < p1; ++j) {
      S11(i,j) = Nmean(i,j) / Ncum[i];
    }
    for (int j{0}; j < p2; ++j) {
      S12(i,j) = Nmean(i,j + p1) / Ncum[i];
    }
  }
  for (int i{0}; i < p2; ++i) {
    for (int j{0}; j < p2; ++j) {
      S22(i,j) = Nmean(i + p1,j + p1) / Ncum[i + p1];
    }
  }
}


//' Loglikelihood for discrete phase-type
//' 
//' @param alpha Initial probabilities.
//' @param S Sub-transition matrix.
//' @param obs The observations.
//' @param weight The weights of the observations.
//' 
// [[Rcpp::export]]
double logLikelihoodDPH(arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight) {
  arma::mat e; 
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = e - (S * e);
  
  double max_val{max(obs)};
  
  std::vector<arma::mat> vect = vector_of_powers(S, max_val);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  for (int k{0}; k < obs.size(); ++k) {
    aux_mat = alpha.t() * vect[obs[k] - 1] * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood for discrete phase-type MoE
//' 
//' @param alpha Initial probabilities.
//' @param S Sub-transition matrix.
//' @param obs The observations.
//' @param weight The weights of the observations.
//' 
// [[Rcpp::export]]
double logLikelihoodDPH_MoE(arma::mat & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight) {
  arma::mat e; 
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = e - (S * e);
  
  double max_val{max(obs)};
  
  std::vector<arma::mat> vect = vector_of_powers(S, max_val);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  for (int k{0}; k < obs.size(); ++k) {
    arma::rowvec alpha_aux(alpha.row(k));
    aux_mat = alpha_aux * vect[obs[k] - 1] * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood for bivariate discrete phase-type
//' 
//' @param alpha Initial probabilities.
//' @param S11 Sub-transition matrix.
//' @param S12 Matrix.
//' @param S22 Sub-transition matrix.
//' @param obs The observations.
//' @param weight The weights of the observations.
//' 
// [[Rcpp::export]]
double logLikelihoodbivDPH(arma::vec & alpha, arma::mat S11, arma::mat S12, arma::mat S22, const Rcpp::NumericMatrix & obs, const Rcpp::NumericVector & weight) {
  arma::mat e; 
  e.ones(S22.n_cols, 1);
  arma::mat exit_vect = e - (S22 * e);
  
  double max_val1{max(obs.column(0))};
  double max_val2{max(obs.column(1))};
  
  std::vector<arma::mat> vect1 = vector_of_powers(S11, max_val1);
  std::vector<arma::mat> vect2 = vector_of_powers(S22, max_val2);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  for (int k{0}; k < obs.nrow(); ++k) {
    aux_mat = alpha.t() * vect1[obs(k, 0) - 1] * S12 * vect2[obs(k, 1) - 1] * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * std::log(density);
  }
  
  return logLh;
}


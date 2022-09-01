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


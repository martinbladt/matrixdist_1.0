# include "EM_LL_UNI.h"
#include "m_exp.h"
# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]

// RIGHT CENSORED 

//' Loglikelihood of matrix-Weibull using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens censored observations.
 //' @param rcweight Weights of the censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMweibull_UNI_inhom(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   
   vector_of_matrices(aux_vect, S, a, m);
   
   arma::mat aux_mat(1,1);
   
   double density{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{pow(obs[k], beta[k])};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density =  aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(beta[k]) + (beta[k] - 1) * std::log(obs[k]));
   }
   // Right censored data
   for (int k{0}; k < rcens.size(); ++k) {
     double x{pow(rcens[k], beta[k])};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * e;
     density = aux_mat(0,0);
     logLh += rcweight[k] * std::log(density);
   }
   
   return logLh;
 }


//' Loglikelihood of matrix-Pareto using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens censored observations.
 //' @param rcweight Weights of the censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMpareto_UNI_inhom(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   
   vector_of_matrices(aux_vect, S, a, m);
   
   arma::mat aux_mat(1,1);
   
   double density{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{std::log(obs[k] / beta[k] + 1)};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) - std::log(obs[k] + beta[k]));
   }
   // Right censored data
   for (int k{0}; k < rcens.size(); ++k) {
     double x{std::log(rcens[k] / beta[k] + 1)};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * e;
     density = aux_mat(0,0);
     logLh += rcweight[k] * std::log(density);
   }
   
   return logLh;
 }


//' Loglikelihood of matrix-lognormal using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens censored observations.
 //' @param rcweight Weights of the censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMlognormal_UNI_inhom(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   
   vector_of_matrices(aux_vect, S, a, m);
   
   arma::mat aux_mat(1,1);
   
   double density{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{pow(std::log(obs[k] + 1), beta[k])};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(beta[k]) + (beta[k] -1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
   }
   // Right censored data
   for (int k{0}; k < rcens.size(); ++k) {
     double x{pow(std::log(rcens[k] + 1), beta[k])};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * e;
     density = aux_mat(0,0);
     logLh += rcweight[k] * std::log(density);
   }
   
   return logLh;
 }


//' Loglikelihood of matrix-loglogistic using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens censored observations.
 //' @param rcweight Weights of the censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMloglogistic_UNI_inhom(double h, arma::vec & alpha, arma::mat & S, arma::mat beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
   if (arma::any(arma::vectorise(beta) < 0)) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   
   vector_of_matrices(aux_vect, S, a, m);
   
   arma::mat aux_mat(1,1);
   
   double density{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{std::log(pow(obs[k] / beta(k,0), beta(k,1)) + 1)};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(beta(k,1)) - std::log(beta(k,0)) + (beta(k,1) - 1) * (std::log(obs[k]) - std::log(beta(k,0))) - std::log(pow(obs[k] / beta(k,0), beta(k,1)) + 1));
   }
   // Right censored data
   for (int k{0}; k < rcens.size(); ++k) {
     double x{std::log(pow(rcens[k] / beta(k,0), beta(k,1)) + 1)};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * e;
     density = aux_mat(0,0);
     logLh += rcweight[k] * std::log(density);
   }
   
   return logLh;
 }


//' Loglikelihood of matrix-Gompertz using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens censored observations.
 //' @param rcweight Weights of the censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMgompertz_UNI_inhom(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   
   vector_of_matrices(aux_vect, S, a, m);
   
   arma::mat aux_mat(1,1);
   
   double density{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{(exp(obs[k] * beta[k]) - 1) / beta[k]};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + obs[k] * beta[k]);
   }
   // Right censored data
   for (int k{0}; k < rcens.size(); ++k) {
     double x{(exp(rcens[k] * beta[k]) - 1) / beta[k]};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * e;
     density = aux_mat(0,0);
     logLh += rcweight[k] * std::log(density);
   }
   
   return logLh;
 }


//' Loglikelihood of matrix-GEV using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens censored observations.
 //' @param rcweight Weights of the censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMgev_UNI_inhom(double h, arma::vec & alpha, arma::mat & S, arma::mat beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
   if (arma::any(arma::vectorise(beta) < 0)) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   
   vector_of_matrices(aux_vect, S, a, m);
   
   arma::mat aux_mat(1,1);
   
   double density{0.0};
   
   double logLh{0.0};
   
   if (arma::all(beta.col(2) == 0)) {
     // Non censored data
     for (int k{0}; k < obs.size(); ++k) {
       double x{exp(-(obs[k] - beta(k,0)) / beta(k,1))};
       
       if (x * a <= 1.0) {
         expm = m_exp_sum(x, m, aux_vect, a);
       }
       else {
         int n{};
         n = std::log(a * x) / std::log(2.0);
         ++n;
         
         expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
         
         pow2_matrix(n, expm);
       }
       aux_mat = alpha.t() * expm * exit_vect;
       density = aux_mat(0,0);
       logLh += weight[k] * (std::log(density) - std::log(beta(k,1)) - (obs[k] - beta(k,0)) / beta(k,1));
     }
     // Right censored data
     for (int k{0}; k < rcens.size(); ++k) {
       double x{exp(-(rcens[k] - beta(k,0)) / beta(k,1))};
       
       if (x * a <= 1.0) {
         expm = m_exp_sum(x, m, aux_vect, a);
       }
       else {
         int n{};
         n = std::log(a * x) / std::log(2.0);
         ++n;
         
         expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
         
         pow2_matrix(n, expm);
       }
       aux_mat = alpha.t() * expm * e;
       density = aux_mat(0,0);
       logLh += rcweight[k] * std::log(density);
     }
   }
   else {
     // Non censored data
     for (int k{0}; k < obs.size(); ++k) {
       double x{pow(1 + (beta(k,2) / beta(k,1)) * (obs[k] - beta(k,0)) , - 1 / beta(k,2))};
       
       if (x * a <= 1.0) {
         expm = m_exp_sum(x, m, aux_vect, a);
       }
       else {
         int n{};
         n = std::log(a * x) / std::log(2.0);
         ++n;
         
         expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
         
         pow2_matrix(n, expm);
       }
       aux_mat = alpha.t() * expm * exit_vect;
       density = aux_mat(0,0);
       logLh += weight[k] * (std::log(density) - std::log(beta(k,1)) - (1 + 1 / beta(k,2)) * std::log(1 + (beta(k,2) / beta(k,1)) * (obs[k] - beta(k,0))));
     }
     // Right censored data
     for (int k{0}; k < rcens.size(); ++k) {
       double x{pow(1 + (beta(k,2) / beta(k,1)) * (rcens[k] - beta(k,0)) , - 1 / beta(k,2))};
       
       if (x * a <= 1.0) {
         expm = m_exp_sum(x, m, aux_vect, a);
       }
       else {
         int n{};
         n = std::log(a * x) / std::log(2.0);
         ++n;
         
         expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
         
         pow2_matrix(n, expm);
       }
       aux_mat = alpha.t() * expm * e;
       density = aux_mat(0,0);
       logLh += rcweight[k] * std::log(density);
     }
   }
   return logLh;
 }


/////////////////////////////////////////////////////////
// Scaled versions of loglikelihoods (for regression)  //
/////////////////////////////////////////////////////////

//' Loglikelihood of PI with matrix-Weibull using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Censored observations.
 //' @param rcweight Weights of the censored observations.
 //' @param scale1 Scale for observations.
 //' @param scale2 Scale for censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMweibull_UNIs_inhom(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   
   vector_of_matrices(aux_vect, S, a, m);
   
   arma::mat aux_mat(1,1);
   
   double density{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{scale1[k] * pow(obs[k], beta[k])};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta[k]) + (beta[k] -1) * std::log(obs[k]));
   }
   // Right censored data
   for (int k{0}; k < rcens.size(); ++k) {
     double x{scale2[k] * pow(rcens[k], beta[k])};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * e;
     density = aux_mat(0,0);
     logLh += rcweight[k] * std::log(density);
   }
   
   return logLh;
 }


//' Loglikelihood of PI with matrix-Pareto using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Censored observations.
 //' @param rcweight Weights of the censored observations.
 //' @param scale1 Scale for observations.
 //' @param scale2 Scale for censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMpareto_UNIs_inhom(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   
   vector_of_matrices(aux_vect, S, a, m);
   
   arma::mat aux_mat(1,1);
   
   double density{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{scale1[k] * std::log(obs[k] / beta[k] + 1)};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) - std::log(obs[k] + beta[k]));
   }
   // Right censored data
   for (int k{0}; k < rcens.size(); ++k) {
     double x{scale2[k] * std::log(rcens[k] / beta[k] + 1)};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * e;
     density = aux_mat(0,0);
     logLh += rcweight[k] * std::log(density);
   }
   
   return logLh;
 }


//' Loglikelihood of PI with matrix-lognormal using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Censored observations.
 //' @param rcweight Weights of the censored observations.
 //' @param scale1 Scale for observations.
 //' @param scale2 Scale for censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMlognormal_UNIs_inhom(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   
   vector_of_matrices(aux_vect, S, a, m);
   
   arma::mat aux_mat(1,1);
   
   double density{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{scale1[k] * pow(std::log(obs[k] + 1), beta[k])};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta[k]) + (beta[k] - 1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
   }
   // Right censored data
   for (int k{0}; k < rcens.size(); ++k) {
     double x{scale2[k] * pow(std::log(rcens[k] + 1), beta[k])};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * e;
     density = aux_mat(0,0);
     logLh += rcweight[k] * std::log(density);
   }
   
   return logLh;
 }


//' Loglikelihood of PI with matrix-loglogistic using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Censored observations.
 //' @param rcweight Weights of the censored observations.
 //' @param scale1 Scale for observations.
 //' @param scale2 Scale for censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMloglogistic_UNIs_inhom(double h, arma::vec & alpha, arma::mat & S, arma::mat beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if (arma::any(arma::vectorise(beta) < 0)) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   
   vector_of_matrices(aux_vect, S, a, m);
   
   arma::mat aux_mat(1,1);
   
   double density{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{scale1[k] * std::log(pow(obs[k] / beta(k,0), beta(k,1)) + 1)};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta(k,1)) - std::log(beta(k,0)) + (beta(k,1) - 1) * (std::log(obs[k]) - std::log(beta(k,0))) - std::log(pow(obs[k] / beta(k,0), beta(k,1)) + 1));
   }
   // Right censored data
   for (int k{0}; k < rcens.size(); ++k) {
     double x{scale2[k] * std::log(pow(rcens[k] / beta(k,0), beta(k,1)) + 1)};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * e;
     density = aux_mat(0,0);
     logLh += rcweight[k] * std::log(density);
   }
   
   return logLh;
 }


//' Loglikelihood of PI with matrix-Gompertz using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Censored observations.
 //' @param rcweight Weights of the censored observations.
 //' @param scale1 Scale for observations.
 //' @param scale2 Scale for censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMgompertz_UNIs_inhom(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   
   vector_of_matrices(aux_vect, S, a, m);
   
   arma::mat aux_mat(1,1);
   
   double density{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{scale1[k] * (exp(obs[k] * beta[k]) - 1) / beta[k]};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + obs[k] * beta[k]);
   }
   // Right censored data
   for (int k{0}; k < rcens.size(); ++k) {
     double x{scale2[k] * (exp(rcens[k] * beta[k]) - 1) / beta[k]};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * e;
     density = aux_mat(0,0);
     logLh += rcweight[k] * std::log(density);
   }
   
   return logLh;
 }

// INTERVAL CENSORING

////////////////////////////////////////////
// Log-likelihoods
////////////////////////////////////////////

//' Loglikelihood of matrix-Weibull using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Interval censored observations.
 //' @param rcweight Weights of the censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMweibull_UNI_inhom_intCens(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   arma::mat expml(size(S));
   arma::mat expmu(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   std::vector<arma::mat> aux_vectl;
   std::vector<arma::mat> aux_vectu;
   
   vector_of_matrices(aux_vect, S, a, m);
   vector_of_matrices(aux_vectl, S, a, m);
   vector_of_matrices(aux_vectu, S, a, m);
   
   arma::mat aux_mat(1,1);
   arma::mat aux_matl(1,1);
   arma::mat aux_matu(1,1);
   
   double density{0.0};
   double densityl{0.0};
   double densityu{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{pow(obs[k], beta[k])};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(beta[k]) + (beta[k] - 1) * std::log(obs[k]));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{pow(rcens(k,0), beta[k])};
     double xu{pow(rcens(k,1), beta[k])};
     
     if (xl * a <= 1.0) {
       expml = m_exp_sum(xl, m, aux_vectl, a);
     }
     else {
       int n{};
       n = std::log(a * xl) / std::log(2.0);
       ++n;
       
       expml = m_exp_sum(xl / pow(2.0, n), m, aux_vectl, a);
       
       pow2_matrix(n, expml);
     }
     
     if (xu * a <= 1.0) {
       expmu = m_exp_sum(xu, m, aux_vectu, a);
     }
     else {
       int n{};
       n = std::log(a * xu) / std::log(2.0);
       ++n;
       
       expmu = m_exp_sum(xu / pow(2.0, n), m, aux_vectu, a);
       
       pow2_matrix(n, expmu);
     }
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1e-5;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of matrix-Pareto using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Interval censored observations.
 //' @param rcweight Weights of the censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMpareto_UNI_inhom_intCens(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   arma::mat expml(size(S));
   arma::mat expmu(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   std::vector<arma::mat> aux_vectl;
   std::vector<arma::mat> aux_vectu;
   
   vector_of_matrices(aux_vect, S, a, m);
   vector_of_matrices(aux_vectl, S, a, m);
   vector_of_matrices(aux_vectu, S, a, m);
   
   arma::mat aux_mat(1,1);
   arma::mat aux_matl(1,1);
   arma::mat aux_matu(1,1);
   
   double density{0.0};
   double densityl{0.0};
   double densityu{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{std::log(obs[k] / beta[k] + 1)};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) - std::log(obs[k] + beta[k]));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{std::log(rcens(k,0) / beta[k] + 1)};
     double xu{std::log(rcens(k,1) / beta[k] + 1)};
     
     if (xl * a <= 1.0) {
       expml = m_exp_sum(xl, m, aux_vectl, a);
     }
     else {
       int n{};
       n = std::log(a * xl) / std::log(2.0);
       ++n;
       
       expml = m_exp_sum(xl / pow(2.0, n), m, aux_vectl, a);
       
       pow2_matrix(n, expml);
     }
     
     if (xu * a <= 1.0) {
       expmu = m_exp_sum(xu, m, aux_vectu, a);
     }
     else {
       int n{};
       n = std::log(a * xu) / std::log(2.0);
       ++n;
       
       expmu = m_exp_sum(xu / pow(2.0, n), m, aux_vectu, a);
       
       pow2_matrix(n, expmu);
     }
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1e-5;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   return logLh;
 }


//' Loglikelihood of matrix-lognormal using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Interval censored observations.
 //' @param rcweight Weights of the censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMlognormal_UNI_inhom_intCens(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   arma::mat expml(size(S));
   arma::mat expmu(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   std::vector<arma::mat> aux_vectl;
   std::vector<arma::mat> aux_vectu;
   
   vector_of_matrices(aux_vect, S, a, m);
   vector_of_matrices(aux_vectl, S, a, m);
   vector_of_matrices(aux_vectu, S, a, m);
   
   arma::mat aux_mat(1,1);
   arma::mat aux_matl(1,1);
   arma::mat aux_matu(1,1);
   
   double density{0.0};
   double densityl{0.0};
   double densityu{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{pow(std::log(obs[k] + 1), beta[k])};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(beta[k]) + (beta[k] -1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{pow(std::log(rcens(k,0) + 1), beta[k])};
     double xu{pow(std::log(rcens(k,1) + 1), beta[k])};
     
     if (xl * a <= 1.0) {
       expml = m_exp_sum(xl, m, aux_vectl, a);
     }
     else {
       int n{};
       n = std::log(a * xl) / std::log(2.0);
       ++n;
       
       expml = m_exp_sum(xl / pow(2.0, n), m, aux_vectl, a);
       
       pow2_matrix(n, expml);
     }
     
     if (xu * a <= 1.0) {
       expmu = m_exp_sum(xu, m, aux_vectu, a);
     }
     else {
       int n{};
       n = std::log(a * xu) / std::log(2.0);
       ++n;
       
       expmu = m_exp_sum(xu / pow(2.0, n), m, aux_vectu, a);
       
       pow2_matrix(n, expmu);
     }
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1e-5;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of matrix-loglogistic using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Interval censored observations.
 //' @param rcweight Weights of the censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMloglogistic_UNI_inhom_intCens(double h, arma::vec & alpha, arma::mat & S, const arma::mat & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight) {
   if (arma::any(arma::vectorise(beta) < 0)) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   arma::mat expml(size(S));
   arma::mat expmu(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   std::vector<arma::mat> aux_vectl;
   std::vector<arma::mat> aux_vectu;
   
   vector_of_matrices(aux_vect, S, a, m);
   vector_of_matrices(aux_vectl, S, a, m);
   vector_of_matrices(aux_vectu, S, a, m);
   
   arma::mat aux_mat(1,1);
   arma::mat aux_matl(1,1);
   arma::mat aux_matu(1,1);
   
   double density{0.0};
   double densityl{0.0};
   double densityu{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{std::log(pow(obs[k] / beta(k,0), beta(k,1)) + 1)};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(beta(k,1)) - std::log(beta(k,0)) + (beta(k,1) - 1) * (std::log(obs[k]) - std::log(beta(k,0))) - std::log(pow(obs[k] / beta(k,0), beta(k,1)) + 1));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{std::log(pow(rcens(k,0) / beta(k,0), beta(k,1)) + 1)};
     double xu{std::log(pow(rcens(k,1) / beta(k,0), beta(k,1)) + 1)};
     
     if (xl * a <= 1.0) {
       expml = m_exp_sum(xl, m, aux_vectl, a);
     }
     else {
       int n{};
       n = std::log(a * xl) / std::log(2.0);
       ++n;
       
       expml = m_exp_sum(xl / pow(2.0, n), m, aux_vectl, a);
       
       pow2_matrix(n, expml);
     }
     
     if (xu * a <= 1.0) {
       expmu = m_exp_sum(xu, m, aux_vectu, a);
     }
     else {
       int n{};
       n = std::log(a * xu) / std::log(2.0);
       ++n;
       
       expmu = m_exp_sum(xu / pow(2.0, n), m, aux_vectu, a);
       
       pow2_matrix(n, expmu);
     }
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1e-5;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of matrix-Gompertz using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Interval censored observations.
 //' @param rcweight Weights of the censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMgompertz_UNI_inhom_intCens(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   arma::mat expml(size(S));
   arma::mat expmu(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   std::vector<arma::mat> aux_vectl;
   std::vector<arma::mat> aux_vectu;
   
   vector_of_matrices(aux_vect, S, a, m);
   vector_of_matrices(aux_vectl, S, a, m);
   vector_of_matrices(aux_vectu, S, a, m);
   
   arma::mat aux_mat(1,1);
   arma::mat aux_matl(1,1);
   arma::mat aux_matu(1,1);
   
   double density{0.0};
   double densityl{0.0};
   double densityu{0.0};
   
   double logLh{0.0};
   double lowLim = 1e-5;
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{(exp(obs[k] * beta[k]) - 1) / beta[k]};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     if(density < lowLim){ density = 1e-5;}
     logLh += weight[k] * (std::log(density) + obs[k] * beta[k]);
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{(exp(rcens(k,0) * beta[k]) - 1) / beta[k]};
     double xu{(exp(rcens(k,1) * beta[k]) - 1) / beta[k]};
     
     if (xl * a <= 1.0) {
       expml = m_exp_sum(xl, m, aux_vectl, a);
     }
     else {
       int n{};
       n = std::log(a * xl) / std::log(2.0);
       ++n;
       
       expml = m_exp_sum(xl / pow(2.0, n), m, aux_vectl, a);
       
       pow2_matrix(n, expml);
     }
     
     if (xu * a <= 1.0) {
       expmu = m_exp_sum(xu, m, aux_vectu, a);
     }
     else {
       int n{};
       n = std::log(a * xu) / std::log(2.0);
       ++n;
       
       expmu = m_exp_sum(xu / pow(2.0, n), m, aux_vectu, a);
       
       pow2_matrix(n, expmu);
     }
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of matrix-GEV using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Interval censored observations.
 //' @param rcweight Weights of the censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMgev_UNI_inhom_intCens(double h, arma::vec & alpha, arma::mat & S, const arma::mat & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight) {
   if(arma::any(beta.col(1) < 0)) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   arma::mat expml(size(S));
   arma::mat expmu(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   std::vector<arma::mat> aux_vectl;
   std::vector<arma::mat> aux_vectu;
   
   vector_of_matrices(aux_vect, S, a, m);
   vector_of_matrices(aux_vectl, S, a, m);
   vector_of_matrices(aux_vectu, S, a, m);
   
   arma::mat aux_mat(1,1);
   arma::mat aux_matl(1,1);
   arma::mat aux_matu(1,1);
   
   double density{0.0};
   double densityl{0.0};
   double densityu{0.0};
   
   double logLh{0.0};
   double lowLim = 1e-5;
   if (arma::all( beta.col(2) == 0)) {
     // Non censored data
     for (int k{0}; k < obs.size(); ++k) {
       double x{exp(-(obs[k] - beta(k,0)) / beta(k,1))};
       
       if (x * a <= 1.0) {
         expm = m_exp_sum(x, m, aux_vect, a);
       }
       else {
         int n{};
         n = std::log(a * x) / std::log(2.0);
         ++n;
         
         expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
         
         pow2_matrix(n, expm);
       }
       aux_mat = alpha.t() * expm * exit_vect;
       density = aux_mat(0,0);
       logLh += weight[k] * (std::log(density) - std::log(beta(k,1)) - (obs[k] - beta(k,0)) / beta(k,1));
     }
     // Interval censored data
     for (int k{0}; k < rcens.n_rows; ++k) {
       
       double xl{exp(-(rcens(k,0) - beta(k,0)) / beta(k,1))};
       double xu{exp(-(rcens(k,1) - beta(k,0)) / beta(k,1))};
       
       if (xl * a <= 1.0) {
         expml = m_exp_sum(xl, m, aux_vectl, a);
       }
       else {
         int n{};
         n = std::log(a * xl) / std::log(2.0);
         ++n;
         
         expml = m_exp_sum(xl / pow(2.0, n), m, aux_vectl, a);
         
         pow2_matrix(n, expml);
       }
       
       if (xu * a <= 1.0) {
         expmu = m_exp_sum(xu, m, aux_vectu, a);
       }
       else {
         int n{};
         n = std::log(a * xu) / std::log(2.0);
         ++n;
         
         expmu = m_exp_sum(xu / pow(2.0, n), m, aux_vectu, a);
         
         pow2_matrix(n, expmu);
       }
       
       aux_matl = alpha.t() * expml * e;
       aux_matu = alpha.t() * expmu * e;
       densityl = aux_matl(0,0);
       densityu = aux_matu(0,0);
       
       double probInt = densityl- densityu;
       if(probInt<lowLim){
         probInt = lowLim;
       }
       logLh += rcweight[k] * std::log(probInt);
     }
   }
   else {
     // Non censored data
     for (int k{0}; k < obs.size(); ++k) {
       double x{pow(1 + (beta(k,2) / beta(k,1)) * (obs[k] - beta(0,k)) , - 1 / beta(k,2))};
       
       if (x * a <= 1.0) {
         expm = m_exp_sum(x, m, aux_vect, a);
       }
       else {
         int n{};
         n = std::log(a * x) / std::log(2.0);
         ++n;
         
         expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
         
         pow2_matrix(n, expm);
       }
       aux_mat = alpha.t() * expm * exit_vect;
       density = aux_mat(0,0);
       logLh += weight[k] * (std::log(density) - std::log(beta(k,1)) - (1 + 1 / beta(k,2)) * std::log(1 + (beta(k,2) / beta(k,1)) * (obs[k] - beta(k,0))));
     }
     // Interval censored data
     for (int k{0}; k < rcens.n_rows; ++k) {
       
       double xl{pow(1 + (beta(k,2) / beta(k,1)) * (rcens(k,0) - beta(k,0)) , - 1 / beta(k,2))};
       double xu{pow(1 + (beta(k,2) / beta(k,1)) * (rcens(k,1) - beta(k,0)) , - 1 / beta(k,2))};
       
       if (xl * a <= 1.0) {
         expml = m_exp_sum(xl, m, aux_vectl, a);
       }
       else {
         int n{};
         n = std::log(a * xl) / std::log(2.0);
         ++n;
         
         expml = m_exp_sum(xl / pow(2.0, n), m, aux_vectl, a);
         
         pow2_matrix(n, expml);
       }
       
       if (xu * a <= 1.0) {
         expmu = m_exp_sum(xu, m, aux_vectu, a);
       }
       else {
         int n{};
         n = std::log(a * xu) / std::log(2.0);
         ++n;
         
         expmu = m_exp_sum(xu / pow(2.0, n), m, aux_vectu, a);
         
         pow2_matrix(n, expmu);
       }
       
       aux_matl = alpha.t() * expml * e;
       aux_matu = alpha.t() * expmu * e;
       densityl = aux_matl(0,0);
       densityu = aux_matu(0,0);
       
       double probInt = densityl- densityu;
       if(probInt<lowLim){
         probInt = lowLim;
       }
       logLh += rcweight[k] * std::log(probInt);
     }
   }
   return logLh;
 }


/////////////////////////////////////////////////////////
// Scaled versions of loglikelihoods (for regression)  //
/////////////////////////////////////////////////////////

//' Loglikelihood of PI with matrix-Weibull using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Interval censored observations.
 //' @param rcweight Weights of the censored observations.
 //' @param scale1 Scale for observations.
 //' @param scale2 Scale for censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMweibull_UNIs_inhom_intCens(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   arma::mat expml(size(S));
   arma::mat expmu(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   std::vector<arma::mat> aux_vectl;
   std::vector<arma::mat> aux_vectu;
   
   vector_of_matrices(aux_vect, S, a, m);
   vector_of_matrices(aux_vectl, S, a, m);
   vector_of_matrices(aux_vectu, S, a, m);
   
   arma::mat aux_mat(1,1);
   arma::mat aux_matl(1,1);
   arma::mat aux_matu(1,1);
   
   double density{0.0};
   double densityl{0.0};
   double densityu{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{scale1[k] * pow(obs[k], beta[k])};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta[k]) + (beta[k] -1) * std::log(obs[k]));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{scale2[k] * pow(rcens(k,0), beta[k])};
     double xu{scale2[k] * pow(rcens(k,1), beta[k])};
     
     if (xl * a <= 1.0) {
       expml = m_exp_sum(xl, m, aux_vectl, a);
     }
     else {
       int n{};
       n = std::log(a * xl) / std::log(2.0);
       ++n;
       
       expml = m_exp_sum(xl / pow(2.0, n), m, aux_vectl, a);
       
       pow2_matrix(n, expml);
     }
     
     if (xu * a <= 1.0) {
       expmu = m_exp_sum(xu, m, aux_vectu, a);
     }
     else {
       int n{};
       n = std::log(a * xu) / std::log(2.0);
       ++n;
       
       expmu = m_exp_sum(xu / pow(2.0, n), m, aux_vectu, a);
       
       pow2_matrix(n, expmu);
     }
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1e-5;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of PI with matrix-Pareto using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Interval censored observations.
 //' @param rcweight Weights of the censored observations.
 //' @param scale1 Scale for observations.
 //' @param scale2 Scale for censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMpareto_UNIs_inhom_intCens(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   arma::mat expml(size(S));
   arma::mat expmu(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   std::vector<arma::mat> aux_vectl;
   std::vector<arma::mat> aux_vectu;
   
   vector_of_matrices(aux_vect, S, a, m);
   vector_of_matrices(aux_vectl, S, a, m);
   vector_of_matrices(aux_vectu, S, a, m);
   
   arma::mat aux_mat(1,1);
   arma::mat aux_matl(1,1);
   arma::mat aux_matu(1,1);
   
   double density{0.0};
   double densityl{0.0};
   double densityu{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{scale1[k] * std::log(obs[k] / beta[k] + 1)};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) - std::log(obs[k] + beta[k]));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{scale2[k] * std::log(rcens(k,0) / beta[k] + 1)};
     double xu{scale2[k] * std::log(rcens(k,1) / beta[k] + 1)};
     
     if (xl * a <= 1.0) {
       expml = m_exp_sum(xl, m, aux_vectl, a);
     }
     else {
       int n{};
       n = std::log(a * xl) / std::log(2.0);
       ++n;
       
       expml = m_exp_sum(xl / pow(2.0, n), m, aux_vectl, a);
       
       pow2_matrix(n, expml);
     }
     
     if (xu * a <= 1.0) {
       expmu = m_exp_sum(xu, m, aux_vectu, a);
     }
     else {
       int n{};
       n = std::log(a * xu) / std::log(2.0);
       ++n;
       
       expmu = m_exp_sum(xu / pow(2.0, n), m, aux_vectu, a);
       
       pow2_matrix(n, expmu);
     }
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1e-5;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
     
   }
   
   return logLh;
 }


//' Loglikelihood of PI with matrix-lognormal using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Interval censored observations.
 //' @param rcweight Weights of the censored observations.
 //' @param scale1 Scale for observations.
 //' @param scale2 Scale for censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMlognormal_UNIs_inhom_intCens(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   arma::mat expml(size(S));
   arma::mat expmu(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   std::vector<arma::mat> aux_vectl;
   std::vector<arma::mat> aux_vectu;
   
   vector_of_matrices(aux_vect, S, a, m);
   vector_of_matrices(aux_vectl, S, a, m);
   vector_of_matrices(aux_vectu, S, a, m);
   
   arma::mat aux_mat(1,1);
   arma::mat aux_matl(1,1);
   arma::mat aux_matu(1,1);
   
   double density{0.0};
   double densityl{0.0};
   double densityu{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{scale1[k] * pow(std::log(obs[k] + 1), beta[k])};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta[k]) + (beta[k] - 1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{scale2[k] * pow(std::log(rcens(k,0) + 1), beta[k])};
     double xu{scale2[k] * pow(std::log(rcens(k,1) + 1), beta[k])};
     
     if (xl * a <= 1.0) {
       expml = m_exp_sum(xl, m, aux_vectl, a);
     }
     else {
       int n{};
       n = std::log(a * xl) / std::log(2.0);
       ++n;
       
       expml = m_exp_sum(xl / pow(2.0, n), m, aux_vectl, a);
       
       pow2_matrix(n, expml);
     }
     
     if (xu * a <= 1.0) {
       expmu = m_exp_sum(xu, m, aux_vectu, a);
     }
     else {
       int n{};
       n = std::log(a * xu) / std::log(2.0);
       ++n;
       
       expmu = m_exp_sum(xu / pow(2.0, n), m, aux_vectu, a);
       
       pow2_matrix(n, expmu);
     }
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1e-5;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of PI with matrix-loglogistic using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Interval censored observations.
 //' @param rcweight Weights of the censored observations.
 //' @param scale1 Scale for observations.
 //' @param scale2 Scale for censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMloglogistic_UNIs_inhom_intCens(double h, arma::vec & alpha, arma::mat & S, const arma::mat & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if (arma::any(arma::vectorise(beta) < 0)) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   arma::mat expml(size(S));
   arma::mat expmu(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   std::vector<arma::mat> aux_vectl;
   std::vector<arma::mat> aux_vectu;
   
   vector_of_matrices(aux_vect, S, a, m);
   vector_of_matrices(aux_vectl, S, a, m);
   vector_of_matrices(aux_vectu, S, a, m);
   
   arma::mat aux_mat(1,1);
   arma::mat aux_matl(1,1);
   arma::mat aux_matu(1,1);
   
   double density{0.0};
   double densityl{0.0};
   double densityu{0.0};
   
   double logLh{0.0};
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{scale1[k] * std::log(pow(obs[k] / beta(k,0), beta(k,1)) + 1)};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta(k,1)) - std::log(beta(k,0)) + (beta(k,1) - 1) * (std::log(obs[k]) - std::log(beta(k,0))) - std::log(pow(obs[k] / beta(k,0), beta(k,1)) + 1));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{scale2[k] * std::log(pow(rcens(k,0) / beta(k,0), beta(k,1)) + 1)};
     double xu{scale2[k] * std::log(pow(rcens(k,1) / beta(k,0), beta(k,1)) + 1)};
     
     if (xl * a <= 1.0) {
       expml = m_exp_sum(xl, m, aux_vectl, a);
     }
     else {
       int n{};
       n = std::log(a * xl) / std::log(2.0);
       ++n;
       
       expml = m_exp_sum(xl / pow(2.0, n), m, aux_vectl, a);
       
       pow2_matrix(n, expml);
     }
     
     if (xu * a <= 1.0) {
       expmu = m_exp_sum(xu, m, aux_vectu, a);
     }
     else {
       int n{};
       n = std::log(a * xu) / std::log(2.0);
       ++n;
       
       expmu = m_exp_sum(xu / pow(2.0, n), m, aux_vectu, a);
       
       pow2_matrix(n, expmu);
     }
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1e-5;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of PI with matrix-Gompertz using uniformization, when the intensity function is regressed on covariate informatino
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta Parameter of transformation.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Interval censored observations.
 //' @param rcweight Weights of the censored observations.
 //' @param scale1 Scale for observations.
 //' @param scale2 Scale for censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodMgompertz_UNIs_inhom_intCens(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if (Rcpp::is_true(Rcpp::any(beta < 0))) return NA_REAL;
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   arma::mat expml(size(S));
   arma::mat expmu(size(S));
   
   double a = max_diagonal(S * (-1));
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   std::vector<arma::mat> aux_vectl;
   std::vector<arma::mat> aux_vectu;
   
   vector_of_matrices(aux_vect, S, a, m);
   vector_of_matrices(aux_vectl, S, a, m);
   vector_of_matrices(aux_vectu, S, a, m);
   
   arma::mat aux_mat(1,1);
   arma::mat aux_matl(1,1);
   arma::mat aux_matu(1,1);
   
   double density{0.0};
   double densityl{0.0};
   double densityu{0.0};
   
   double logLh{0.0};
   double lowLim = 1e-5;
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{scale1[k] * (exp(obs[k] * beta[k]) - 1) / beta[k]};
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     density = aux_mat(0,0);
     if(density < lowLim){ density = 1e-5;}
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + obs[k] * beta[k]);
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{scale2[k] * (exp(rcens(k,0) * beta[k]) - 1) / beta[k]};
     double xu{scale2[k] * (exp(rcens(k,1) * beta[k]) - 1) / beta[k]};
     
     
     if (xl * a <= 1.0) {
       expml = m_exp_sum(xl, m, aux_vectl, a);
     }
     else {
       int n{};
       n = std::log(a * xl) / std::log(2.0);
       ++n;
       
       expml = m_exp_sum(xl / pow(2.0, n), m, aux_vectl, a);
       
       pow2_matrix(n, expml);
     }
     
     if (xu * a <= 1.0) {
       expmu = m_exp_sum(xu, m, aux_vectu, a);
     }
     else {
       int n{};
       n = std::log(a * xu) / std::log(2.0);
       ++n;
       
       expmu = m_exp_sum(xu / pow(2.0, n), m, aux_vectu, a);
       
       pow2_matrix(n, expmu);
     }
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     
     if(probInt<lowLim){
       probInt = 1e-5;
     }
     
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }

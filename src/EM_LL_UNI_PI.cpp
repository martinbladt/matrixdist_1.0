#include "EM_LL_UNI.h"
#include "m_exp.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Loglikelihood of PI using uniformization, with customizable observation transformation.
 //'
 //' Loglikelihood for a sample, automatically handling right-censored or interval-censored data based on the structure of `rcens`.
 //'
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param beta1 Parameter of transformation (uncensored observations).
 //' @param beta2 Parameter of transformation (censored observations).
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Censored observations (either numeric vector for right censoring or matrix for interval censoring).
 //' @param rcweight Weights of the censored observations.
 //' @param scale1 Scale for observations.
 //' @param scale2 Scale for censored observations.
 //' @param gfun_name Name of the transformation function to be applied to observations (e.g., "scale", "log", "sqrt").
 //' 
 // [[Rcpp::export]]
 double logLikelihood_UNIs_PI(double h, arma::vec & alpha,
                              arma::mat & S,
                              SEXP beta1,
                              SEXP beta2,
                              const Rcpp::NumericVector & obs,
                              const Rcpp::NumericVector & weight,
                              SEXP rcens, 
                              const Rcpp::NumericVector & rcweight,
                              const Rcpp::NumericVector & scale1,
                              const Rcpp::NumericVector & scale2,
                              const std::string & gfun_name) {
   
   // Declare beta as either NumericVector or NumericMatrix depending on beta1's type
   bool is_beta1_vector = Rcpp::is<Rcpp::NumericVector>(beta1);
   bool is_beta2_vector = Rcpp::is<Rcpp::NumericVector>(beta2);
   
   bool is_beta1_matrix = Rcpp::is<Rcpp::NumericMatrix>(beta1);
   bool is_beta2_matrix = Rcpp::is<Rcpp::NumericMatrix>(beta2);
   
   if (!is_beta1_vector && !is_beta1_matrix) {
     Rcpp::stop("beta1 must be either a NumericVector or NumericMatrix");
   }
   if (!is_beta2_vector && !is_beta2_matrix) {
     Rcpp::stop("beta1 must be either a NumericVector or NumericMatrix");
   }
   
   
   // Cast beta1 into the appropriate type
   NumericVector beta1_vec;
   NumericMatrix beta1_mat;
   NumericVector beta2_vec;
   NumericMatrix beta2_mat;
   
   if (is_beta1_vector) {
     beta1_vec = as<NumericVector>(beta1);
     if (Rcpp::is_true(Rcpp::any(beta1_vec < 0))) {
       return NA_REAL;
     }
   } else if (is_beta1_matrix) {
     beta1_mat = as<NumericMatrix>(beta1);
     if (Rcpp::is_true(Rcpp::any(beta1_mat < 0))) {
       return NA_REAL;
     }
   }
   
   if (is_beta2_vector) {
     beta2_vec = as<NumericVector>(beta2);
     if (Rcpp::is_true(Rcpp::any(beta2_vec < 0))) {
       return NA_REAL;
     }
   } else if (is_beta2_matrix) {
     beta2_mat = as<NumericMatrix>(beta2);
     if (Rcpp::is_true(Rcpp::any(beta2_mat < 0))) {
       return NA_REAL;
     }
   }
   
   arma::mat e;
   e.ones(S.n_cols, 1);
   arma::mat exit_vect = (S * (-1)) * e;
   
   arma::mat expm(size(S));
   arma::mat expml(size(S)), expmu(size(S));
   
   double a = max_diagonal(S * (-1));
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   vector_of_matrices(aux_vect, S, a, m);
   
   arma::mat aux_mat(1, 1);
   
   double density{0.0};
   double logLh{0.0};
   
   // Non-censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x;
     if (gfun_name == "weibull") {
       x = scale1[k] * pow(obs[k], beta1_vec[k]);
       
     } else if (gfun_name == "pareto") {
       x = scale1[k] * std::log(obs[k] / beta1_vec[k] + 1);
       
     } else if (gfun_name == "lognormal") {
       x = scale1[k] * pow(std::log(obs[k] + 1), beta1_vec[k]);
       
     } else if (gfun_name == "loglogistic") {
       x = scale1[k] * std::log(pow(obs[k] / beta1_mat(k,0), beta1_mat(k,1)) + 1);
       
     }  else if (gfun_name == "gompertz") {
       x = scale1[k] * (exp(obs[k] * beta1_vec[k]) - 1) / beta1_vec[k];
       
     }else {
       Rcpp::stop("Unknown gfun_name: " + gfun_name);
     }
     
     if (x * a <= 1.0) {
       expm = m_exp_sum(x, m, aux_vect, a);
     } else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       pow2_matrix(n, expm);
     }
     aux_mat = alpha.t() * expm * exit_vect;
     double density = aux_mat(0, 0);
     
     // Determine the logarithm of the density function
     if (gfun_name == "weibull") {
       logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta1_vec[k]) + (beta1_vec[k] - 1) * std::log(obs[k]));
       
     } else if (gfun_name == "pareto") {
       logLh += weight[k] * (std::log(density) + std::log(scale1[k]) - std::log(obs[k] + beta1_vec[k]));
       
     } else if (gfun_name == "lognormal") {
       logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta1_vec[k]) + (beta1_vec[k] - 1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
       
     } else if (gfun_name == "loglogistic") {
       logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta1_mat(k,1)) - std::log(beta1_mat(k,0)) + (beta1_mat(k,1) - 1) * (std::log(obs[k]) - std::log(beta1_mat(k,0))) - std::log(pow(obs[k] / beta1_mat(k,0), beta1_mat(k,1)) + 1));
       
     }  else if (gfun_name == "gompertz") {
       logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + obs[k] * beta1_vec[k]);
     }
   }
   
   // Check if rcens is a NumericVector (right-censored) or NumericMatrix (interval-censored)
   if (Rcpp::is<Rcpp::NumericVector>(rcens)) {
     Rcpp::NumericVector rcens_vec(rcens);
     
     // Right-censored data
     for (int k{0}; k < rcens_vec.size(); ++k) {
       double x;
       if (gfun_name == "weibull") {
         x = scale2[k] * pow(rcens_vec[k], beta2_vec[k]);
         
       } else if (gfun_name == "pareto") {
         x = scale2[k] * std::log(rcens_vec[k] / beta2_vec[k] + 1);
         
       } else if (gfun_name == "lognormal") {
         x = scale2[k] * pow(std::log(rcens_vec[k] + 1), beta2_vec[k]);
         
       } else if (gfun_name == "loglogistic") {
         x = scale2[k] * std::log(pow(rcens_vec[k] / beta2_mat(k,0), beta2_mat(k,1)) + 1);
         
       }  else if (gfun_name == "gompertz") {
         x = scale2[k] * (exp(rcens_vec[k] * beta2_vec[k]) - 1) / beta2_vec[k];
         
       }
       
       if (x * a <= 1.0) {
         expm = m_exp_sum(x, m, aux_vect, a);
       } else {
         int n{};
         n = std::log(a * x) / std::log(2.0);
         ++n;
         
         expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
         pow2_matrix(n, expm);
       }
       aux_mat = alpha.t() * expm * e;
       double density = aux_mat(0, 0);
       logLh += rcweight[k] * std::log(density);
     }
     
   } else if (Rcpp::is<Rcpp::NumericMatrix>(rcens)) {
     Rcpp::NumericMatrix rcens_mat(rcens);
     
     // Interval-censored data
     std::vector<arma::mat> aux_vectl, aux_vectu;
     vector_of_matrices(aux_vectl, S, a, m);
     vector_of_matrices(aux_vectu, S, a, m);
     
     for (int k{0}; k < rcens_mat.nrow(); ++k) {
       double xl, xu;
       if (gfun_name == "weibull") {
         double xl{pow(rcens_mat(k,0), beta2_vec[k])};
         double xu{pow(rcens_mat(k,1), beta2_vec[k])};
         
       } else if (gfun_name == "pareto") {
         double xl{std::log(rcens_mat(k,0) / beta2_vec[k] + 1)};
         double xu{std::log(rcens_mat(k,1) / beta2_vec[k] + 1)};
         
       } else if (gfun_name == "lognormal") {
         double xl{pow(std::log(rcens_mat(k,0) + 1), beta2_vec[k])};
         double xu{pow(std::log(rcens_mat(k,1) + 1), beta2_vec[k])};
         
       } else if (gfun_name == "loglogistic") {
         double xl{std::log(pow(rcens_mat(k,0) / beta2_mat(k,0), beta2_mat(k,1)) + 1)};
         double xu{std::log(pow(rcens_mat(k,1) / beta2_mat(k,0), beta2_mat(k,1)) + 1)};
         
       } else if (gfun_name == "gompertz") {
         double xl{(exp(rcens_mat(k,0) * beta2_vec[k]) - 1) / beta2_vec[k]};
         double xu{(exp(rcens_mat(k,1) * beta2_vec[k]) - 1) / beta2_vec[k]};
         
       }
       
       // Lower bound of interval
       if (xl * a <= 1.0) {
         expml = m_exp_sum(xl, m, aux_vectl, a);
       } else {
         int n{};
         n = std::log(a * xl) / std::log(2.0);
         ++n;
         
         expml = m_exp_sum(xl / pow(2.0, n), m, aux_vectl, a);
         pow2_matrix(n, expml);
       }
       
       // Upper bound of interval
       if (xu * a <= 1.0) {
         expmu = m_exp_sum(xu, m, aux_vectu, a);
       } else {
         int n{};
         n = std::log(a * xu) / std::log(2.0);
         ++n;
         
         expmu = m_exp_sum(xu / pow(2.0, n), m, aux_vectu, a);
         pow2_matrix(n, expmu);
       }
       
       arma::mat aux_matl = alpha.t() * expml * e;
       arma::mat aux_matu = alpha.t() * expmu * e;
       double densityl = aux_matl(0, 0);
       double densityu = aux_matu(0, 0);
       
       double probInt = densityl - densityu;
       if (probInt < 1e-5) probInt = 1e-5;
       
       logLh += rcweight[k] * std::log(probInt);
     }
   } else {
     Rcpp::stop("rcens must be either a numeric vector (right-censored) or numeric matrix (interval-censored).");
   }
   
   return logLh;
 }
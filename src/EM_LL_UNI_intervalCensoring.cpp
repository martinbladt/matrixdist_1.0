# include "EM_LL_UNI.h"
#include "m_exp.h"
# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]

//' EM for phase-type using uniformization for matrix exponential
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param obs The observations.
 //' @param weight The weights for the observations.
 //' @param rcens Matrix of censored observations' bounds. Each row corresponds to a given observation.
 //' @param rcweight The weights for the censored observations.
 //' 
 // [[Rcpp::export]]
 void EMstep_UNI_intervalCensoring(double h,
                                   arma::vec & alpha,
                                   arma::mat & S,
                                   const Rcpp::NumericVector & obs,
                                   const Rcpp::NumericVector & weight,
                                   const arma::mat & rcens,
                                   const Rcpp::NumericVector & rcweight) {
   unsigned p{S.n_rows};
   
   // column unit vector
   arma::mat e;
   e.ones(S.n_cols, 1); 
   // pxp matrix of zeros
   arma::mat zeroMat;
   zeroMat.zeros(p,p);
   // column exit rate vector
   arma::mat exit_vect = (S * (-1)) * e;
   // Identity matrix
   arma::mat id;
   id.ones(p,p);
   
   // sufficient statistics
   arma::mat Bmean = arma::zeros(p,1);
   arma::mat Zmean = arma::zeros(p,1);
   arma::mat Nmean = arma::zeros(p,p + 1);
   
   // row vector of {alpha*exp(S*x)}
   arma::mat avector(1,p);
   // col vector of {exp(S*x)*s}
   arma::mat bvector(p,1);
   // Van Loan matrix of {\int_0^x \exp{(S(x-u))} e \alpha \exp{(Su)} du}
   arma::mat cmatrix(p,p);
   // matrix exponential {\exp{Sx}}
   arma::mat aux_exp(p,p);
   
   // same as above for interval lower bound
   arma::mat bvectorl(p,1);
   arma::mat cmatrixl(p,p);
   arma::mat cmatrixlint(p,p);
   arma::mat aux_expl(p,p);
   
   // same as above for interval upper bound
   arma::mat bvectoru(p,1);
   arma::mat cmatrixu(p,p);
   arma::mat cmatrixuint(p,p);
   arma::mat aux_expu(p,p);
   
   // PH density
   arma::mat aux_mat(1,1);
   // prob of surviving the lower bound
   arma::mat aux_matl(1,1);
   // prob of surviving the upper bound
   arma::mat aux_matu(1,1);
   
   // \int_l^u \exp{Su} du
   arma::mat intluMat(p,p);
   // \alpha(\int_l^u \exp{Su} du)e_k
   arma::mat intluVec(1,p);
   // {\int_a^b \exp{(S(a-u))} e \alpha \exp{(Su)} du}
   arma::mat intluVanLoan(p,p);
   
   // block matrices for Van Loan matrix exponentials
   arma::mat J(2 * p,2 * p);
   arma::mat Jl(2 * p,2 * p);
   arma::mat Ju(2 * p,2 * p);
   arma::mat Jlintegral(2 * p,2 * p);
   arma::mat Juintegral(2 * p,2 * p);
   
   // matrix for Van Loan {e\alpha}
   arma::mat s_prod_alpha(p,p);
   s_prod_alpha = exit_vect * alpha.t();
   
   J = matrix_vanloan(S, S, s_prod_alpha);
   
   double a = max_diagonal(J * (-1));
   double aInt;
   
   int m{find_n(h, 1)};
   
   std::vector<arma::mat> aux_vect;
   std::vector<arma::mat> aux_vectl;
   std::vector<arma::mat> aux_vectu;
   
   std::vector<arma::mat> aux_vectlintegral;
   std::vector<arma::mat> aux_vectuintegral;
   
   vector_of_matrices(aux_vect, J, a, m);
   
   double sum_weights{0.0};
   double density{0.0};
   double densityl{0.0};
   double densityu{0.0};
   
   // E-step
   // Uncensored data
   for (int k{0}; k < obs.size(); ++k) {
     sum_weights += weight[k];
     
     double x{obs[k]};
     
     if (x * a <= 1.0) {
       J = m_exp_sum(x, m, aux_vect, a);
     }
     else {
       int n{};
       n = std::log(a * x) / std::log(2.0);
       ++n;
       
       J = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
       
       pow2_matrix(n, J);
     }
     
     for (int i{0}; i < p; ++i) {
       for (int j{0}; j < p; ++j) {
         aux_exp(i,j) = J(i,j);
         cmatrix(i,j) = J(i,j + p);
       }
     }
     
     avector = alpha.t() * aux_exp;
     bvector = aux_exp * exit_vect;
     aux_mat = alpha.t() * bvector;
     density = aux_mat(0,0);
     
     // E-step
     for (int i{0}; i < p; ++i) {
       Bmean(i,0) += alpha[i] * bvector(i,0) * weight[k] / density;
       Nmean(i,p) += avector(0,i) * exit_vect(i,0) * weight[k] / density;
       Zmean(i,0) += cmatrix(i,i) * weight[k] / density;
       for (int j{0}; j < p; ++j) {
         Nmean(i,j) += S(i,j) * cmatrix(j,i) * weight[k] / density;
       }
     }
   }
   
   // Interval-Censored Data
   double sum_censored{0.0};
   if (rcens.n_rows> 0) {
     
     s_prod_alpha = e * alpha.t();
     
     Jl = matrix_vanloan(S, S, s_prod_alpha);
     Ju = matrix_vanloan(S, S, s_prod_alpha);
     
     Jlintegral = matrix_vanloan(S, zeroMat, id);
     Juintegral = matrix_vanloan(S, zeroMat, id);
     
     a = max_diagonal(Jl * (-1));
     
     aInt = max_diagonal(Jlintegral * (-1));
     
     aux_vect.clear();
     vector_of_matrices(aux_vectl, Jl, a, m);
     vector_of_matrices(aux_vectu, Ju, a, m);
     
     vector_of_matrices(aux_vectlintegral, Jlintegral, a, m);
     vector_of_matrices(aux_vectuintegral, Juintegral, a, m);
   }
   for (int k{0}; k < rcens.n_rows; ++k) {
     sum_censored += rcweight[k];
     
     // Lower bound
     double xl{rcens(k,0)};
     if (xl * a <= 1.0) {
       Jl = m_exp_sum(xl, m, aux_vectl, a);
       Jlintegral = m_exp_sum(xl, m, aux_vectlintegral, a);
       
     } else {
       int nl{};
       int nlInt{};
       
       nl = std::log(a * xl) / std::log(2.0);
       ++nl;
       
       nlInt = std::log(a * xl) / std::log(2.0);
       ++nlInt;
       
       Jl = m_exp_sum(xl / pow(2.0, nl), m, aux_vectl, a);
       Jlintegral = m_exp_sum(xl / pow(2.0, nlInt), m, aux_vectlintegral, a);
       
       pow2_matrix(nl, Jl);
       pow2_matrix(nlInt, Jlintegral);
     }
     
     // Upper bound
     double xu{rcens(k,1)};
     if (xu * a <= 1.0) {
       Ju = m_exp_sum(xu, m, aux_vectu, a);
       Juintegral = m_exp_sum(xu, m, aux_vectuintegral, a); 
       
     } else {
       int nu{};
       int nuInt{};
       
       nu = std::log(a * xu) / std::log(2.0);
       ++nu;
       
       nuInt = std::log(a * xu) / std::log(2.0);
       ++nuInt;
       
       Ju = m_exp_sum(xu / pow(2.0, nu), m, aux_vectu, a);
       Juintegral = m_exp_sum(xu / pow(2.0, nuInt), m, aux_vectuintegral, a);
       
       pow2_matrix(nu, Ju);
       pow2_matrix(nuInt, Juintegral);
     }
     
     for (int i{0}; i < p; ++i) {
       for (int j{0}; j < p; ++j) {
         // Lower bound
         aux_expl(i,j) = Jl(i,j);
         cmatrixl(i,j) = Jl(i,j + p);
         cmatrixlint(i,j) = Jlintegral(i,j + p);
         
         // Upper bound
         aux_expu(i,j) = Ju(i,j);
         cmatrixu(i,j) = Ju(i,j + p);
         cmatrixuint(i,j) = Juintegral(i,j + p);
       }
     }
     
     // Lower bound
     bvectorl = aux_expl * e;
     aux_matl = alpha.t() * bvectorl;
     densityl = aux_matl(0,0);
     
     // Upper bound
     bvectoru = aux_expu * e;
     aux_matu = alpha.t() * bvectoru;
     densityu = aux_matu(0,0);
     
     intluMat = cmatrixuint-cmatrixlint;
     intluVec = alpha.t()*intluMat;
     intluVanLoan = cmatrixu-cmatrixl;
     
     // E-step
     for (int i{0}; i < p; ++i) {
       Bmean(i,0) += alpha[i] * (bvectorl(i,0)-bvectoru(i,0)) * rcweight[k] / (densityl-densityu);
       Zmean(i,0) += (intluVec(0,i)-intluVanLoan(i,i))* rcweight[k] / (densityl-densityu);
       Nmean(i,p) += exit_vect(i,0) * intluVec(0,i) * rcweight[k] / (densityl-densityu);
       for (int j{0}; j < p; ++j) {
         Nmean(i,j) += S(i,j) * (intluVec(0,i)-intluVanLoan(j,i)) * rcweight[k] / (densityl-densityu);
       }
     }
   }
   
   // M step
   for (int i{0}; i < p; ++i) {
     alpha[i] = Bmean(i,0) / (sum_weights + sum_censored);
     if (alpha[i] < 0) {
       alpha[i] = 0;
     }
     exit_vect(i,0) = Nmean(i,p) / Zmean(i,0);
     if (exit_vect(i,0) < 0) {
       exit_vect(i,0) = 0;
     }
     S(i,i) = -exit_vect(i,0);
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

//' Loglikelihood of phase-type using uniformization
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Interval censored observations.
 //' @param rcweight Weights of the censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodPH_UNI_intervalCensoring(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight) {
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
     double x{obs[k]};
     
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
     logLh += weight[k] * std::log(density);
   }
   //Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     double xl{rcens(k,0)};
     double xu{rcens(k,1)};
     
     expml = matrix_exponential(S*xl);
     expmu = matrix_exponential(S*xu);
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl-densityu;
     double lowLim = 1/100000;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of matrix-Weibull using uniformization
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
 double logLikelihoodMweibull_UNI_intervalCensoring(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight) {
   if(beta < 0) return NA_REAL;
   
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
     double x{pow(obs[k], beta)};
     
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
     logLh += weight[k] * (std::log(density) + std::log(beta) + (beta - 1) * std::log(obs[k]));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{pow(rcens(k,0), beta)};
     double xu{pow(rcens(k,1), beta)};
     
     expml = matrix_exponential(S*xl);
     expmu = matrix_exponential(S*xu);
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1/100000;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of matrix-Pareto using uniformization
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
 double logLikelihoodMpareto_UNI_intervalCensoring(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight) {
   if(beta < 0) return NA_REAL;
   
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
     double x{std::log(obs[k] / beta + 1)};
     
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
     logLh += weight[k] * (std::log(density) - std::log(obs[k] + beta));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{std::log(rcens(k,0) / beta + 1)};
     double xu{std::log(rcens(k,1) / beta + 1)};
     
     expml = matrix_exponential(S*xl);
     expmu = matrix_exponential(S*xu);
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     double probInt = densityl- densityu;
     double lowLim = 1/100000;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   return logLh;
 }


//' Loglikelihood of matrix-lognormal using uniformization
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
 double logLikelihoodMlognormal_UNI_intervalCensoring(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight) {
   if(beta < 0) return NA_REAL;
   
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
     double x{pow(std::log(obs[k] + 1), beta)};
     
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
     logLh += weight[k] * (std::log(density) + std::log(beta) + (beta -1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{pow(std::log(rcens(k,0) + 1), beta)};
     double xu{pow(std::log(rcens(k,1) + 1), beta)};
     
     expml = matrix_exponential(S*xl);
     expmu = matrix_exponential(S*xu);
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1/100000;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of matrix-loglogistic using uniformization
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
 double logLikelihoodMloglogistic_UNI_intervalCensoring(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight) {
   if(beta[0] < 0 || beta[1] < 0) return NA_REAL;
   
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
     double x{std::log(pow(obs[k] / beta[0], beta[1]) + 1)};
     
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
     logLh += weight[k] * (std::log(density) + std::log(beta[1]) - std::log(beta[0]) + (beta[1] - 1) * (std::log(obs[k]) - std::log(beta[0])) - std::log(pow(obs[k] / beta[0], beta[1]) + 1));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{std::log(pow(rcens(k,0) / beta[0], beta[1]) + 1)};
     double xu{std::log(pow(rcens(k,1) / beta[0], beta[1]) + 1)};
     
     expml = matrix_exponential(S*xl);
     expmu = matrix_exponential(S*xu);
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1/100000;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of matrix-Gompertz using uniformization
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
 double logLikelihoodMgompertz_UNI_intervalCensoring(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight) {
   if(beta < 0) return NA_REAL;
   
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
   double lowLim = 1/100000;
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{(exp(obs[k] * beta) - 1) / beta};
     
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
     if(density < lowLim){ density = 1/100000;}
     logLh += weight[k] * (std::log(density) + obs[k] * beta);
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{(exp(rcens(k,0) * beta) - 1) / beta};
     double xu{(exp(rcens(k,1) * beta) - 1) / beta};
     
     expml = matrix_exponential(S*xl);
     expmu = matrix_exponential(S*xu);
     
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


//' Loglikelihood of matrix-GEV using uniformization
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
 double logLikelihoodMgev_UNI_intervalCensoring(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight) {
   if(beta[1] < 0) return NA_REAL;
   
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
   double lowLim = 1/100000;
   if (beta[2] == 0) {
     // Non censored data
     for (int k{0}; k < obs.size(); ++k) {
       double x{exp(-(obs[k] - beta[0]) / beta[1])};
       
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
       logLh += weight[k] * (std::log(density) - std::log(beta[1]) - (obs[k] - beta[0]) / beta[1]);
     }
     // Interval censored data
     for (int k{0}; k < rcens.n_rows; ++k) {
       
       double xl{exp(-(rcens(k,0) - beta[0]) / beta[1])};
       double xu{exp(-(rcens(k,1) - beta[0]) / beta[1])};
       
       expml = matrix_exponential(S*xl);
       expmu = matrix_exponential(S*xu);
       
       aux_matl = alpha.t() * expml * e;
       aux_matu = alpha.t() * expmu * e;
       densityl = aux_matl(0,0);
       densityu = aux_matu(0,0);
       
       double probInt = densityl- densityu;
       double lowLim = 1/100000;
       if(probInt<lowLim){
         probInt = lowLim;
       }
       logLh += rcweight[k] * std::log(probInt);
     }
   }
   else {
     // Non censored data
     for (int k{0}; k < obs.size(); ++k) {
       double x{pow(1 + (beta[2] / beta[1]) * (obs[k] - beta[0]) , - 1 / beta[2])};
       
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
       logLh += weight[k] * (std::log(density) - std::log(beta[1]) - (1 + 1 / beta[2]) * std::log(1 + (beta[2] / beta[1]) * (obs[k] - beta[0])));
     }
     // Interval censored data
     for (int k{0}; k < rcens.n_rows; ++k) {
       
       double xl{pow(1 + (beta[2] / beta[1]) * (rcens(k,0) - beta[0]) , - 1 / beta[2])};
       double xu{pow(1 + (beta[2] / beta[1]) * (rcens(k,1) - beta[0]) , - 1 / beta[2])};
       
       expml = matrix_exponential(S*xl);
       expmu = matrix_exponential(S*xu);
       
       aux_matl = alpha.t() * expml * e;
       aux_matu = alpha.t() * expmu * e;
       densityl = aux_matl(0,0);
       densityu = aux_matu(0,0);
       
       double probInt = densityl- densityu;
       double lowLim = 1/100000;
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

//' Loglikelihood of PI with phase-type using uniformization
 //' 
 //' Loglikelihood for a sample.
 //' 
 //' @param h Positive parameter.
 //' @param alpha Initial probabilities.
 //' @param S Sub-intensity matrix.
 //' @param obs The observations.
 //' @param weight Weights of the observations.
 //' @param rcens Interval censored observations.
 //' @param rcweight Weights of the censored observations.
 //' @param scale1 Scale for observations.
 //' @param scale2 Scale for censored observations.
 //' 
 // [[Rcpp::export]]
 double logLikelihoodPH_UNIs_intervalCensoring(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
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
     double x{scale1[k] * obs[k]};
     
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
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{scale2[k] * rcens(k,0)};
     double xu{scale2[k] * rcens(k,1)};
     
     expml = matrix_exponential(S*xl);
     expmu = matrix_exponential(S*xu);
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1/100000;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of PI with matrix-Weibull using uniformization
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
 double logLikelihoodMweibull_UNIs_intervalCensoring(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if(beta < 0) return NA_REAL;
   
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
     double x{scale1[k] * pow(obs[k], beta)};
     
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
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta) + (beta -1) * std::log(obs[k]));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{scale2[k] * pow(rcens(k,0), beta)};
     double xu{scale2[k] * pow(rcens(k,1), beta)};
     
     expml = matrix_exponential(S*xl);
     expmu = matrix_exponential(S*xu);
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1/100000;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of PI with matrix-Pareto using uniformization
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
 double logLikelihoodMpareto_UNIs_intervalCensoring(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if(beta < 0) return NA_REAL;
   
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
     double x{scale1[k] * std::log(obs[k] / beta + 1)};
     
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
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) - std::log(obs[k] + beta));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{scale2[k] * std::log(rcens(k,0) / beta + 1)};
     double xu{scale2[k] * std::log(rcens(k,1) / beta + 1)};
     
     expml = matrix_exponential(S*xl);
     expmu = matrix_exponential(S*xu);
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1/100000;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
     
   }
   
   return logLh;
 }


//' Loglikelihood of PI with matrix-lognormal using uniformization
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
 double logLikelihoodMlognormal_UNIs_intervalCensoring(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if(beta < 0) return NA_REAL;
   
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
     double x{scale1[k] * pow(std::log(obs[k] + 1), beta)};
     
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
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta) + (beta - 1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{scale2[k] * pow(std::log(rcens(k,0) + 1), beta)};
     double xu{scale2[k] * pow(std::log(rcens(k,1) + 1), beta)};
     
     expml = matrix_exponential(S*xl);
     expmu = matrix_exponential(S*xu);
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1/100000;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of PI with matrix-loglogistic using uniformization
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
 double logLikelihoodMloglogistic_UNIs_intervalCensoring(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if(beta[0] < 0 || beta[1] < 0) return NA_REAL;
   
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
     double x{scale1[k] * std::log(pow(obs[k] / beta[0], beta[1]) + 1)};
     
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
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta[1]) - std::log(beta[0]) + (beta[1] - 1) * (std::log(obs[k]) - std::log(beta[0])) - std::log(pow(obs[k] / beta[0], beta[1]) + 1));
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{scale2[k] * std::log(pow(rcens(k,0) / beta[0], beta[1]) + 1)};
     double xu{scale2[k] * std::log(pow(rcens(k,1) / beta[0], beta[1]) + 1)};
     
     expml = matrix_exponential(S*xl);
     expmu = matrix_exponential(S*xu);
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     double lowLim = 1/100000;
     if(probInt<lowLim){
       probInt = lowLim;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }


//' Loglikelihood of PI with matrix-Gompertz using Uniformization
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
 double logLikelihoodMgompertz_UNIs_intervalCensoring(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const arma::mat & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
   if(beta < 0) return NA_REAL;
   
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
   double lowLim = 1/100000;
   
   // Non censored data
   for (int k{0}; k < obs.size(); ++k) {
     double x{scale1[k] * (exp(obs[k] * beta) - 1) / beta};
     
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
     if(density < lowLim){ density = 1/100000;}
     logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + obs[k] * beta);
   }
   // Interval censored data
   for (int k{0}; k < rcens.n_rows; ++k) {
     
     double xl{scale2[k] * (exp(rcens(k,0) * beta) - 1) / beta};
     double xu{scale2[k] * (exp(rcens(k,1) * beta) - 1) / beta};
     
     expml = matrix_exponential(S*xl);
     expmu = matrix_exponential(S*xu);
     
     aux_matl = alpha.t() * expml * e;
     aux_matu = alpha.t() * expmu * e;
     densityl = aux_matl(0,0);
     densityu = aux_matu(0,0);
     
     double probInt = densityl- densityu;
     
     if(probInt<lowLim){
       probInt = 1/100000;
     }
     logLh += rcweight[k] * std::log(probInt);
   }
   
   return logLh;
 }

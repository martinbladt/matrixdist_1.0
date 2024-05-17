# include <RcppArmadillo.h>
# include "m_exp.h"

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
#endif


// [[ Rcpp :: depends ( RcppArmadillo )]]



arma::mat kron_sum2(arma::mat A, arma::mat B) {
  arma::mat I1;
  I1.eye(size(A));
  arma::mat I2;
  I2.eye(size(B));
  arma::mat K = kron(A, I2) + kron(I1, B);
  return K;
}


//' Bivariate phase-type joint density of the common shock type
//'
//' @param x Matrix of values.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-intensity matrix.
//' @param P Matrix.
//' @param Q1 Sub-intensity matrix.
//' @param Q2 Sub-intensity matrix.
//' @return Joint density at \code{x}.
//' @export
//' 
// [[Rcpp::export]]
Rcpp::NumericVector csph_density_par(Rcpp::NumericMatrix x, arma::vec alpha, arma::mat S, arma::mat P, arma::mat Q1, arma::mat Q2) {
   unsigned p1{S.n_rows};
   unsigned p2{Q1.n_rows};
   long n{x.nrow()};
   
   
   Rcpp::NumericVector density(n);
   
   arma::mat e1;
   e1.ones(Q1.n_cols, 1);
   arma::mat exit_vect1 = (Q1 * (-1)) * e1;
   
   arma::mat e2;
   e2.ones(Q2.n_cols, 1);
   arma::mat exit_vect2 = (Q2 * (-1)) * e2;
   
   arma::mat exit_vect_prod = kron(exit_vect1, exit_vect2);
   arma::mat Q1pQ2 = kron_sum2(Q1, Q2);
   
   
   
   arma::mat aux_mat(1,1);
   
  #pragma omp parallel
  {
  #pragma omp for
   for (int k{0}; k < n; ++k) {
     double m = std::min(x(k,0), x(k,1));
     arma::mat B1 = matrix_exponential(Q1 * (x(k,0) - m));
     arma::mat B2 = matrix_exponential(Q2 * (x(k,1) - m));
     arma::mat B1tB2 = kron(B1, B2);
     arma::mat b_prod_alpha = B1tB2 * exit_vect_prod * alpha.t();
     arma::mat cmatrix(p2 * p2, p1);
     
     for (int i{0}; i < p2; ++i) {
       arma::mat ei = arma::zeros(1,p2);
       ei[i] = 1;
       arma::mat eitei = kron(ei, ei);
       arma::mat J = matrix_exponential(matrix_vanloan(Q1pQ2, S, b_prod_alpha) * m);
       for (int l{0}; l < p2 * p2; ++l) {
         for (int j{0}; j < p1; ++j) {
           cmatrix(l,j) = J(l,j + p2 * p2);
         }
       }
       aux_mat = eitei * cmatrix * P * ei.t();
       density[k] += aux_mat(0,0);
     }
     
   }
   
  }
   return density;
 }

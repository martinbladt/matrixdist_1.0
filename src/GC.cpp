# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]

// [[Rcpp::export]]
Rcpp::NumericVector gcdensity(const Rcpp::NumericVector & x, arma::vec & alpha, arma::mat & S) {
  
  long p{S.n_rows};
  
  Rcpp::NumericVector density(x.size());
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::vec t = vectorise(S * e);
  
  arma::vec lambda = diagvec(S *(-1));
  
  arma::vec mu = 1+ t/lambda;
  
  arma::vec pr = lambda % mu;
  
  Rcpp::NumericVector b(x.size());
  
  for(int j = 0; j < p; ++j){
    if(alpha[j] != 0){
      for(int k = j; k < p; ++k){
        double fact{1.0};
        if(k > j) {
          fact = prod(pr.subvec(j, k - 1));
        }
        double a = fact * lambda[k] * (1 - mu[k]);
        b = b * 0;
        for(int m = j; m <= k; ++m){
          arma::vec factor = lambda.subvec(j,k) - lambda[m];
          fact = prod(factor.elem( find(factor != 0) ));
          b = b + exp(-lambda[m] * x) / fact;
        }
        density = density + alpha[j] * a * b;
      }
    }
  }
  
  return density;
}
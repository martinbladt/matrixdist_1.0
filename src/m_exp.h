# include <RcppArmadillo.h>

double inf_norm(arma::mat A);

arma::mat matrix_VanLoan(arma::mat A1, arma::mat A2, arma::mat B1);

double max_diagonal(const arma::mat & A);

arma::mat matrix_exponential(arma::mat A);
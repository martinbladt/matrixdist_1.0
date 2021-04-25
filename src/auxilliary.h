# include <RcppArmadillo.h>

double LInf_normArma(arma::mat A);

arma::mat mexponentialArma(arma::mat Ainput);

arma::mat matrix_VanLoanArma(arma::mat A1, arma::mat A2, arma::mat B1);

double matrixMaxDiagonal_arma(const arma::mat & A);
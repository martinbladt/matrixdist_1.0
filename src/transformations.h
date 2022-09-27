#ifndef TVR   // if x.h hasn't been included yet...
#define TVR   //   #define this so the compiler knows it has been included

#include <RcppArmadillo.h>

int n_pos (const arma::vec R);

int n_null (const arma::vec R);

arma::vec plus_states (const arma::vec R);

arma::vec null_states (const arma::vec R);

arma::mat Q_pos_pos (const arma::vec R, const arma::mat Qtilda);

arma::mat Q_null_null (const arma::vec R, const arma::mat Qtilda);

arma::mat Q_pos_null (const arma::vec R, const arma::mat Qtilda);

arma::mat Q_null_pos (const arma::vec R, const arma::mat Qtilda);

arma::vec q_pos (const arma::vec R, const arma::mat Qtilda);

arma::vec q_null (const arma::vec R, const arma::mat Qtilda);

arma::mat new_trans_mat (const arma::vec R, const arma::mat Qtilda);

arma::vec new_trans_exit (const arma::vec R, const arma::mat Qtilda);

arma::rowvec pi_pos (const arma::vec R, const arma::vec alpha);

arma::rowvec pi_null (const arma::vec R, const arma::vec alpha);

arma::vec new_pi(arma::vec R, arma::mat Qtilda, arma::vec alpha);

arma::vec new_exit_vec (const arma::vec R, const arma::mat Qtilda, const arma::mat S);

arma::mat new_subint_mat (const arma::vec R, const arma::mat Qtilda, const arma::mat S);

Rcpp::List transf_via_rew(arma::mat R,arma::mat Qtilda, arma::vec alpha, arma::mat S );
#endif

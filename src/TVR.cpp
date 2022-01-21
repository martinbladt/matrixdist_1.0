#include <RcppArmadillo.h>
#include "TVR.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

// Intermediary functions, to improve readability
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

//' Find how many states have positive reward
//'
//' @param R reward vector
//'
//' @return The number of states with positive rewards
//'
// [[Rcpp::export]]
int n_pos(arma::vec R) {
  int d{0};
  
  for (int i{0}; i < R.size(); ++i) {
    if (R(i) > 0){
      ++d;
    }
  }
  return d;
}

//' Find how many states have null reward
//'
//' @param R reward vector
//'
//' @return The number of states with null rewards
//'
// [[Rcpp::export]]
int n_null(arma::vec R) {
  int l{0};
  
  for (int i{0}; i < R.size(); ++i) {
    if (R(i) == 0) {
      ++l;
    }
  }
  return l;
}

//' Find which states have positive reward
//'
//' @param R reward vector
//'
//' @return A vector with the states (number) that are associated with positive rewards
//'
// [[Rcpp::export]]
arma::vec plus_states(arma::vec R) {
  int d = n_pos(R);
  int j{0};
  
  arma::vec plusState(d);
  
  for (int i{0}; i < R.size(); ++i) {
    if (R(i) > 0) {
      plusState(j) = i + 1;
      ++j;
    }
  }
  return plusState;
}

//' Find which states have null reward
//'
//' @param R reward vector
//'
//' @return A vector with the states (number) that are associated with null rewards
//'
// [[Rcpp::export]]
arma::vec null_states(arma::vec R) {
  int l = n_null(R);
  int j{0};
  
  arma::vec zeroState(l);
  
  for (int i=0; i < R.size(); ++i) {
    if (R(i) == 0) {
      zeroState(j) = i + 1;
      ++j;
    }
  }
  return zeroState;
}

//' Obtain Q++
//'
//' @param R reward vector
//' @param Qtilda transition-matrix of the Embedded MC
//'
//' @return The matrix Q++ where we have transition probabilities from/to states associated with positive rewards
//'
// [[Rcpp::export]]
arma::mat Q_pos_pos(arma::vec R, arma::mat Qtilda) {
  int n{0};
  int m{0};
  int d = n_pos(R);
  
  arma::mat Qpp(d,d);
  arma::vec plusState = plus_states(R);
  arma::vec copy = plusState;
  
  for (int i:plusState) {
    m = 0;
    for (int j:copy) {
      Qpp(n,m) = Qtilda(i - 1, j - 1);
      ++m;
    }
    ++n;
  }
  return(Qpp);
}

//' Obtain Q00
//'
//' @param R reward vector
//' @param Qtilda transition-matrix of the Embedded MC
//'
//' @return the matrix Q00 where we have transition probabilities from/to states associated with null rewards
//'
// [[Rcpp::export]]
arma::mat Q_null_null (arma::vec R, arma::mat Qtilda) {
  int n{0};
  int m{0};
  int l = n_null(R);
  
  arma::mat Qzz(l,l);
  arma::vec zeroState = null_states(R);
  arma::vec copy = zeroState;
  
  for (int i:zeroState) {
    m = 0;
    for (int j:copy) {
      Qzz(n,m) = Qtilda(i - 1,j - 1);
      ++m;
    }
    ++n;
  }
  return(Qzz);
}

//' Obtain Qtilda+0
//'
//'
//' @param R reward vector
//' @param Qtilda transition-matrix of the Embedded MC
//'
//' @return the matrix Qtilda+0 where we have transition probabilities from states associated with positive rewards to ones associated with null rewards
//'
// [[Rcpp::export]]
arma::mat Q_pos_null(arma::vec R, arma::mat Qtilda) {
  int n{0};
  int m{0};
  int l = n_null(R);
  int d = n_pos(R);
  
  arma::mat Qpz(d,l);
  arma::vec zeroState = null_states(R);
  arma::vec plusState = plus_states(R);
  
  for (int i:plusState) {
    m = 0;
    for (int j:zeroState) {
      Qpz(n,m) = Qtilda(i - 1,j - 1);
      ++m;
    }
    ++n;
  }
  return Qpz;
}

//' Obtain Q0+
//'
//'
//' @param R reward vector
//' @param Qtilda transition-matrix of the Embedded MC
//'
//' @return the matrix Q0+ where we have transition probabilities from states associated with null rewards  to ones associated with positive rewards
//'
// [[Rcpp::export]]
arma::mat Q_null_pos(arma::vec R, arma::mat Qtilda) {
  int n{0};
  int m{0};
  int l = n_null(R);
  int d = n_pos(R);
  
  arma::mat Qzp(l,d);
  arma::vec zeroState = null_states(R);
  arma::vec plusState = plus_states(R);
  
  for (int i:zeroState) {
    m = 0;
    for (int j:plusState) {
      Qzp(n,m) = Qtilda(i - 1,j - 1);
      ++m;
    }
    ++n;
  }
  return Qzp;
}

//' Obtain q+
//'
//' @param R reward vector
//' @param Qtilda transition-matrix of the Embedded MC
//'
//' @return Creates the vector with transition probabilities from states associated with positive rewards to the absorption state
//'
// [[Rcpp::export]]
arma::vec q_pos(arma::vec R, arma::mat Qtilda){
  unsigned p{Qtilda.n_rows};
  int d = n_pos(R);
  
  arma::vec q(p);
  arma::vec qp(d);
  arma::vec plusState = plus_states(R);
  
  q = Qtilda.col(p - 1);
  
  int j{0};
  
  for (int i:plusState) {
    qp(j) = q(i - 1);
    ++j;
  }
  
  return qp;
}

//' Obtain q0
//'
//' @param R reward vector
//' @param Qtilda transition-matrix of the Embedded MC
//'
//' @return Creates the vector with transition probabilities from states associated with null rewards to the absorption state
//'
// [[Rcpp::export]]
arma::vec q_null(arma::vec R, arma::mat Qtilda) {
  unsigned p{Qtilda.n_rows};
  int l = n_null(R);
  
  arma::vec q(p);
  arma::vec qz(l);
  arma::vec zeroState = null_states(R);
  
  q = Qtilda.col(p - 1);
  
  int j{0};
  
  for (int i:zeroState) {
    qz(j) = q(i - 1);
    ++j;
  }
  
  return qz;
}

//' Obtain new sub-transition matrix for the new embedded MC
//'
//' @param R reward vector
//' @param Qtilda transition-matrix of the Embedded MC
//'
//' @return The sub-transition matrix of the new embedded Markov Chain
//'
// [[Rcpp::export]]
arma::mat new_trans_mat(arma::vec R, arma::mat Qtilda) {
  int l = n_null(R);
  
  arma::mat Qpp = Q_pos_pos(R,Qtilda);
  arma::mat Qzz = Q_null_null(R,Qtilda);
  arma::mat Qpz = Q_pos_null(R,Qtilda);
  arma::mat Qzp = Q_null_pos(R,Qtilda);
  arma::mat inter = inv(arma::eye(l,l) - Qzz);
  
  arma::mat P = Qpp + Qpz * inter * Qzp;
  
  return P;
}

//' Obtain new sub-transition exit vector for the new embedded MC
//'
//' @param R reward vector
//' @param Qtilda transition-matrix of the Embedded MC
//'
//' @return The exit rates for the new embedded Markov Chain
//'
// [[Rcpp::export]]
arma::vec new_trans_exit(arma::vec R, arma::mat Qtilda) {
  int d = n_pos(R);
  
  arma::mat P = new_trans_mat(R,Qtilda);
  arma::mat e(d,1); 
  e.fill(1);
  
  arma::mat p = e - P * e;
  
  return p;
}

//' Get pi+
//'
//' @param R reward vector
//' @param alpha initial distribution vector of the original MJP
//'
//' @return The initial distribution vector for states associated with positive rewards
//'
// [[Rcpp::export]]
arma::rowvec pi_pos(arma::vec R,arma::vec alpha){
  int d = n_pos(R);
  
  arma::rowvec piP(d);
  arma::vec plusState = plus_states(R);
  
  int j{0};
  
  for (int i:plusState) {
    piP(j) = alpha(i - 1);
    ++j;
  }
  
  return piP;
}

//' Get pi0
//'
//' @param R reward vector
//' @param alpha initial distribution vector of the original MJP
//'
//' @return The initial distribution vector for states associated with null rewards
//'
// [[Rcpp::export]]
arma::rowvec pi_null(arma::vec R,arma::vec alpha) {
  int l = n_null(R);
  
  arma::rowvec piZ(l);
  arma::vec zeroState = null_states(R);
  
  int j{0};
  
  for (int i:zeroState) {
    piZ(j) = alpha(i-1);
    ++j;
  }
  
  return piZ;
}

//' Obtain new sub-transition matrix for the new embedded MC
//'
//' @param R reward vector
//' @param Qtilda transition-matrix of the Embedded MC
//' @param alpha the original initial distribution vector
//'
//'@return The initial distribution vector for the new Markov Jump Process
//'
// [[Rcpp::export]]
arma::vec new_pi(arma::vec R, arma::mat Qtilda, arma::vec alpha) {
  int l = n_null(R);
  
  arma::rowvec piP = pi_pos(R,alpha);
  arma::rowvec piZ = pi_null(R,alpha);
  arma::mat Qzz = Q_null_null(R,Qtilda);
  arma::mat Qzp = Q_null_pos(R,Qtilda);
  
  arma::mat inter = inv(arma::eye(l,l)-Qzz);
  
  arma::rowvec alphaNew(l);
  
  alphaNew= piP + piZ * inter * Qzp;
  
  return alphaNew.t();
  
}

//' Obtain new exit rate vector for the new embedded MC
//'
//' @param R reward vector
//' @param Qtilda transition-matrix of the Embedded MC
//' @param S original sub-intensity matrix
//'
//' @return The exit rate vector for the new Markov Jump Process
//'
// [[Rcpp::export]]
arma::vec new_exit_vec(arma::vec R, arma::mat Qtilda, arma::mat S) {
  int d = n_pos(R);
  int index{0};
  
  arma::vec plusState = plus_states(R);
  
  arma::vec p = new_trans_exit(R,Qtilda);
  arma::vec tstar(d);
  
  for (int i{0}; i < d; ++i) {
    index = plusState(i) - 1;
    tstar(i) = -S(index,index) / R(index) * p(i);
  }
  
  return tstar;
}

//' Obtain new sub-intensity matrix for the new MJP
//'
//' @param R reward vector
//' @param Qtilda transition-matrix of the Embedded MC
//' @param S original sub-intensity matrix
//'
//' @return The sub-intensity matrix for the new Markov Jump Process
//'
// [[Rcpp::export]]
arma::mat new_subint_mat(arma::vec R, arma::mat Qtilda, arma::mat S) {
  int d = n_pos(R);
  int index{0};
  
  arma::vec plusState = plus_states(R);
  arma::vec zeroState = null_states(R);
  
  arma::mat P = new_trans_mat(R,Qtilda);
  arma::vec tstar = new_exit_vec(R,Qtilda,S);
  
  arma::mat Tstar(d,d);
  
  for (int i{0}; i < d; ++i) {
    for (int j{0}; j < d; ++j) {
      if (j != i) {
        index = plusState(i) - 1;
        Tstar(i,j) = -S(index,index) / R(index) * P(i,j);
      }
    }
  }
  
  double rest{0};
  
  for (int k{0}; k < d; ++k) {
    rest = sum(Tstar.row(k));
    Tstar(k,k) = -(rest + tstar(k));
  }
  
  return Tstar;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//The TVR function that will be called to perform the transformations.

//' Performs the TVR
//'
//' @param R reward vector
//' @param Qtilda transition-matrix of the Embedded MC
//' @param alpha initial distribution vector
//' @param S original sub-intensity matrix
//'
//' @return A list of transformed PH distributions
//'
// [[Rcpp::export]]
Rcpp::List transf_via_rew(arma::mat R,arma::mat Qtilda, arma::vec alpha, arma::mat S ) {
  int k = R.n_cols;
  
  Rcpp::List Marginal;
  
  for (int i{0}; i < k; ++i) {
    arma::vec reward;
    
    arma::vec alpha1;
    arma::mat S1;
    
    reward = R.col(i);
    
    alpha1 = new_pi(reward, Qtilda, alpha);
    if (reward.is_zero() == FALSE) {
      S1 = new_subint_mat(reward, Qtilda, S);
    } else{ 
      S1.zeros(1,1);
    }
    
    Rcpp::List x = Rcpp::List::create(Rcpp::Named("alpha") = alpha1.t(), Rcpp::_["S"] = S1);
    Marginal.push_back(x);
  }
  
  return Marginal;
}

#include <RcppArmadillo.h>
#include "TVR.h"
#include "EM_LL_UNI.h"
#include "Simulation.h"
#include "m_exp.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
using namespace arma;

//' Random multivariate phase-type via rewards
//'
//' Generates samples of size \code{n} from a multivariate phase-type distribution with parameters \code{alpha} and \code{S}.
//' @param n Sample size for each marginal.
//' @param alpha Initial probabilities.
//' @param S sub-intensity matrix.
//' @param R reward matrix.
//'
//' @return The simulated samples, each column corresponds to a marginal.
//'
// [[Rcpp::export]]
arma::mat rMPHstar(int n, arma::vec alpha, arma::mat S, arma::mat R) {
  
  int d{R.n_cols}; //number of marginals
  
  List Marginal=transf_via_rew(R, embedded_mc(S), alpha, S);
  List marginalJ;
  
  arma::mat mph_sample(n,d);
  
  for(int m{0}; m<d ; ++m){
    marginalJ=Marginal(m);
    mph_sample.col(m)=Rcpp::as<arma::vec>(rphasetype(n,marginalJ[0],marginalJ[1]));
  }
  
  return (mph_sample);
}

//' Random reward matrix
//'
//' Generates a random reward matrix for a multivariate phase-type distribution with p states and d marginals
//'
//' @param p  the number of transient states in the sub-intensity matrix
//' @param d  the number of marginals
//'
//' @return The random reward matrix
//'
// [[Rcpp::export]]
arma::mat random_reward(long int p, long int d){
  
  arma::mat R(p,d);
  arma::vec line(d);
  arma::vec tot(d);
  
  for(long int i=0; i<p;++i){
    
    line=arma::randu<arma::vec>(d);
    tot.fill(arma::sum(line));
    line=line/tot;
    
    R.row(i)=line.t();
    
  }
  return(R);
  
}

//' Transform a reward matrix with very small rewards to avoid numerical problems
//'
//' @param R Reward matrix
//' @param tol Lower bound considered for a reward
//'
//' @return A reward matrix that does not cause issues with uniformization
//'
// [[Rcpp::export]]
void rew_sanity_check (arma::mat & R,double tol){
  int p{R.n_rows};
  int d{R.n_cols};
  
  for(int i{0}; i<p; ++i){
    for(int j{0}; j<d; ++j){
      
      if(R(i,j)<tol){
        R(i,j)=0;}
    }
  }
  
  arma::vec e(d); e.ones();
  arma::vec check=R*e;
  
  // if a reward is transformed its value is redistributed to other states with positive rewards,
  // redistribution is based on the weight of remaining rewards.
  
  for(int i{0}; i<check.size(); ++i){
    
    if(check(i)!=1){
      
      arma::vec rew=plus_states(conv_to<vec>::from(R.row(i)));
      
      double miss{1-check(i)};
      
      for(int k{0};k<rew.size();++k){
        int j{static_cast<int>(rew(k))-1};
        double w{R(i,j)/check(i)};
        R(i,j)+=miss*w;
      }
    }
    
  }
}

////////////////////////////////////////////////////////////////
//  Multivariate EM (algoritm 2 in paper) via Uniformization  //
////////////////////////////////////////////////////////////////


//' Marginal conditional expectations
//'
//' @param rew Column of the reward matrix corresponding to its marginal.
//' @param pos Vector that indicates which state is associated to a positive reward.
//' @param N Uniformization parameter.
//' @param obs Marginal observations.
//' @param weight Marginal weights.
//' @param rcens Marginal right-censored values.
//' @param rcweight Marginal weights for rc values.
//' @param alpha Marginal initial distribution vector.
//' @param S Marginal sub-intensity matrix.
//' 
//' @return A vector with the expected time spent in each state by the marginal, conditional on the observations.
//'
//' @export
// [[Rcpp::export]]
arma::vec marginal_expectation (arma::vec & rew, arma::vec & pos,int N, arma::vec & alpha, arma::mat & S,arma::vec & obs, arma::vec & weight, arma::vec & rcens, arma::vec & rcweight ){
  
  arma::vec Zmeanj(rew.size()); Zmeanj.zeros();
  
  int p{S.n_cols};  //the dimension of the mth marginal (can be smaller than p if some rewards are null)
  
  double density{0.0};
  
  //E-step for marginal distributions
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat bvector(p,1);
  arma::mat cmatrix(p,p);
  arma::mat aux_exp(p,p);
  
  arma::mat aux_mat(1,1);
  
  arma::mat J(2 * p,2 * p);
  arma::mat tProductPi(p,p);
  tProductPi = t * alpha.t();
  
  J = matrix_vanloan(S, S, tProductPi);
  
  double a = max_diagonal(J * (-1));
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, J, a, N);
  
  //  Unccensored data
  if(obs.size()>0){
    for (int k{0}; k < obs.size(); ++k) {
      
      double x{obs[k]};
      
      if (x * a <= 1.0) {
        J = m_exp_sum(x, N, aux_vect, a);
      }
      else {
        int n{};
        n = std::log(a * x) / std::log(2.0);
        ++n;
        
        J = m_exp_sum(x / pow(2.0, n), N, aux_vect, a);
        
        pow2_matrix(n, J);
      }
      
      for (int i{0}; i < p; ++i) {
        for (int j{0}; j < p; ++j) {
          aux_exp(i,j) = J(i,j);
          cmatrix(i,j) = J(i,j + p);
        }
      }
      if(aux_exp.is_zero()==TRUE){aux_exp.print("exp(Sx):"); Rcout<< "exp(Sx):"<< endl; stop("Issue with uniformization for marg uncensored data-> exp(Sx) is a matrix of zeros");}
      if(aux_exp.has_nan()==TRUE){aux_exp.print("exp(Sx):"); Rcout<< "exp(Sx):"<< endl; stop("At least one NaN element in exp(Sx) in marg un data");}
      if(cmatrix.has_nan()==TRUE){cmatrix.print("G(x;alpha,S):"); Rcout<< "G(x;alpha,S):" << endl; stop("At least one NaN element in G(x;alpha,S) in marg un data");}
      
      bvector = aux_exp * t;
      aux_mat = alpha.t() * bvector;
      density = aux_mat(0,0);
      
      //E-step for expected time spent in state k by marginal j
      int j{0};
      for (int i{0}; i < p; ++i) {
        j=pos(i)-1;
        Zmeanj(j) += cmatrix(i,i) * weight[k] / density;
      }
    }
  }
  
  
  //  Right-Censored Data for jth marginal
  if (rcens.size()>0) {
    tProductPi = e * alpha.t();
    J = matrix_vanloan(S, S, tProductPi);
    aux_vect.clear();
    vector_of_matrices(aux_vect, J, a, N);
    
    
    for (int k{0}; k < rcens.size(); ++k) {
      
      double x{rcens[k]};
      
      if (x * a <= 1.0) {
        J = m_exp_sum(x, N, aux_vect, a);
      }
      else {
        int n{};
        n = std::log(a * x) / std::log(2.0);
        ++n;
        
        J = m_exp_sum(x / pow(2.0, n), N, aux_vect, a);
        
        pow2_matrix(n, J);
      }
      
      for (int i{0}; i < p; ++i) {
        for (int j{0}; j < p; ++j) {
          aux_exp(i,j) = J(i,j);
          cmatrix(i,j) = J(i,j + p);
        }
      }
      if(aux_exp.is_zero()==TRUE){aux_exp.print("exp(Sx):"); Rcout<< "exp(Sx):"<< endl; stop("Issue with uniformization for marg rc data-> exp(Sx) is a matrix of zeros");}
      if(aux_exp.has_nan()==TRUE){aux_exp.print("exp(Sx):"); Rcout<< "exp(Sx):"<< endl; stop("At least one NaN element in exp(Sx) for marg rc data");}
      if(cmatrix.has_nan()==TRUE){cmatrix.print("G(x;alpha,S):"); Rcout<< "G(x;alpha,S):" << endl; stop("At least one NaN element in G(x;alpha,S) for marg rc");}
      
      bvector = aux_exp * e;
      aux_mat = alpha.t() * bvector;
      density = aux_mat(0,0);
      
      //E-step
      int j{0};
      for (int i{0}; i < p; ++i) {
        j=pos(i)-1;
        Zmeanj(j) += cmatrix(i,i) * rcweight[k] / density;
      }
    }
  }
  
  return(Zmeanj);
}

//' EM step using Uniformization for MPHstar class
//'
//' @param h positive parameter for precision of uniformization method.
//' @param Rtol The smallest value that a reward can take.
//' @param alpha Vector of initial probabilities of the originating distribution.
//' @param S The sub-intensity matrix of the originating distribution.
//' @param R The reward matrix.
//' @param mph_obs The list of summed, marginal observations (uncensored and right censored) with associated weights.
//'
//' @export
// [[Rcpp::export]]
void MPHstar_EMstep_UNI(double h, double Rtol, arma::vec & alpha, arma::mat & S,arma::mat & R ,const Rcpp::List & mph_obs) {
  
  rew_sanity_check(R, Rtol); // to avoid numerical instability we consider a reward 0 if it is smaller than Rtol
  
  int p{S.n_cols}; //dimension of original PH distribution
  int d{R.n_cols}; //number of marginal distributions
  
  Rcpp::List preobs=mph_obs[0]; //List of observations and weights for summed data
  arma::mat un_sum=preobs["un"];
  arma::mat rc_sum=preobs["rc"];
  
  arma::vec un_obsSum=un_sum.col(0); // vector of unique sum of observations (uncensored)
  arma::vec unWeightSum=un_sum.col(1); // vector of weights for unique observations (uncensored)
  
  arma::vec rcensSum=rc_sum.col(0); // vector of unique right censoring values
  arma::vec rcweightSum=rc_sum.col(1); // vector of weights for censoring values
  
  arma::mat Qtilda=embedded_mc(S); // probability matrix of the embedded MC
  
  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat ed; ed.ones(d,1);
  arma::mat t = (S * (-1)) * e;
  
  arma::mat Bmean = arma::zeros(p,1);
  arma::mat Zmean = arma::zeros(p,1);
  arma::mat Nmean = arma::zeros(p,p + 1);
  
  arma::mat ZmeanTot = arma::zeros(p,1);
  arma::mat Zmeanj = arma::zeros(p,d); // for E[Z^j_k | Y^j=y^j], rows for states and columns for marginal
  
  // Uniformization matrices
  arma::mat avector(1,p);
  arma::mat bvector(p,1);
  arma::mat cmatrix(p,p);
  arma::mat aux_exp(p,p);
  
  arma::mat aux_mat(1,1);
  
  arma::mat J(2 * p,2 * p);
  arma::mat tProductPi(p,p);
  tProductPi = t * alpha.t();
  
  J = matrix_vanloan(S, S, tProductPi); // Van Loan (1978)  get the block matrix J=(S,tProductPi;0,S)
  
  double a = max_diagonal(J * (-1));
  
  int N{find_n(h, 1)}; //  Find n such that P(N > n) = h with N Poisson(1) distributed
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, J, a, N); // Computes elements P^n / n! until the value size N, P=I +J*1/a
  // with P=I + J *(1/ a), we have aux_vect=(I P P^2/(2!) ... P^N/(N!) )
  
  
  double SumOfWeights{0.0};
  double density{0.0};
  
  //E-step for sum of marginal data
  //  Unccensored data
  for (int k{0}; k < un_obsSum.size(); ++k) {
    SumOfWeights += unWeightSum[k];
    
    double x{un_obsSum[k]};
    
    if (x * a <= 1.0) {
      J = m_exp_sum(x, N, aux_vect, a); // J= I + P*(ax) +(P(ax))^2/2!+...+(P(ax))^N/N!
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      J = m_exp_sum(x / pow(2.0, n), N, aux_vect, a); // J= I + P*(x/n^2a)+...+ (P(ax/n^2))^N/N!
      
      pow2_matrix(n, J); // J^n -> now J=(exp(Sx),G(x;pi,S);0,exp(Sx)) 2px2p block matrix, according to Van Loan(1978)
    }
    
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);  // aux_exp=exp(Sx)
        cmatrix(i,j) = J(i,j + p); // cmatrix=G(x;pi,S)
      }
    }
    if(aux_exp.is_zero()==TRUE){aux_exp.print("exp(Sx):"); Rcout<< "exp(Sx):"<< endl; stop("Issue with uniformization for summed unc. data-> exp(Sx) is a matrix of zeros");}
    if(aux_exp.has_nan()==TRUE){aux_exp.print("exp(Sx):"); Rcout<< "exp(Sx):"<< endl; stop("At least one NaN element in exp(Sx) for summed unc. data");}
    if(cmatrix.has_nan()==TRUE){cmatrix.print("G(x;alpha,S):"); Rcout<< "G(x;alpha,S):" << endl; stop("At least one NaN element in G(x;alpha,S) for summed unc.");}
    
    
    avector = alpha.t() * aux_exp; // avector= pi%*%exp(Sx)
    bvector = aux_exp * t; // bvector=exp(Sx)%*%t
    aux_mat = alpha.t() * bvector; //aux_mat=pi%*%exp(Sx)%*%t ->density of original PH
    density = aux_mat(0,0); // density of original PH
    
    //E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha(i) * bvector(i,0) * unWeightSum[k] / density;
      Nmean(i,p) += avector(0,i) * t(i,0) *unWeightSum[k] / density;
      Zmean(i,0) += cmatrix(i,i) * unWeightSum[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += S(i,j) * cmatrix(j,i) * unWeightSum[k] / density;
      }
    }
  }
  
  //  Right-Censored Data
  double SumOfCensored{0.0};
  
  if (rcensSum.size() > 0) {
    tProductPi = e * alpha.t(); // for rc data we have e%*%pi instead of t%*%pi
    J = matrix_vanloan(S, S, tProductPi);
    aux_vect.clear();
    vector_of_matrices(aux_vect, J, a, N);
  }
  for (int k{0}; k < rcensSum.size(); ++k) {
    SumOfCensored += rcweightSum(k);
    
    double x{rcensSum(k)}; //observation is now the right censoring value
    
    //Uniformization same as for Uncensored observations
    if (x * a <= 1.0) {
      J = m_exp_sum(x, N, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      J = m_exp_sum(x / pow(2.0, n), N, aux_vect, a);
      
      pow2_matrix(n, J);
    }
    
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }
    if(aux_exp.is_zero()==TRUE){aux_exp.print("exp(Sx):"); Rcout<< "exp(Sx):"<< endl; stop("Issue with uniformization for summed rc data-> exp(Sx) is a matrix of zeros");}
    if(aux_exp.has_nan()==TRUE){aux_exp.print("exp(Sx):"); Rcout<< "exp(Sx):"<< endl; stop("At least one NaN element in exp(Sx) for summed rc data");}
    if(cmatrix.has_nan()==TRUE){cmatrix.print("G(x;alpha,S):"); Rcout<< "G(x;alpha,S):" << endl; stop("At least one NaN element in G(x;alpha,S) for summed rc data");}
    
    
    bvector = aux_exp * e; //bvector=exp(Sx)%*%e
    aux_mat = alpha.t() * bvector; // aux_mat=pi%*%exp(Sx)%*%e -> P(tau>x)
    density = aux_mat(0,0); //
    
    //E-step
    //No expectation for N_k since we consider that absorption did not occur yet
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha(i) * bvector(i,0) * rcweightSum(k) / density;
      Zmean(i,0) += cmatrix(i,i) * rcweightSum(k) / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += S(i,j) * cmatrix(j,i) * rcweightSum(k) / density;
      }
    }
  }
  
  
  // E-step for marginals
  
  for(int m{0};m<d;++m){
    arma::vec rew=R.col(m);
    arma::vec pos=plus_states(rew);
    
    int pj{pos.size()};
    
    if(pj==0){continue;
    }else if(pj==1){
      arma::mat alphaj=conv_to<mat>::from(new_pi(rew, Qtilda, alpha)); // initial distribution vector of mth marginal
      arma::mat Sj=new_subint_mat(rew, Qtilda, S); // sub-intensity matrix of mth marginal
      
      for(int i{0}; i<rew.size();++i){
        if(rew(i)>0){Zmeanj(i,m)=alphaj(0,0)*1/Sj(0,0);}
      }
      
    }else{
      
      arma::vec alphaj=new_pi(rew, Qtilda, alpha); // initial distribution vector of mth marginal
      arma::mat Sj=new_subint_mat(rew, Qtilda, S); // sub-intensity matrix of mth marginal
      
      
      preobs=mph_obs[m+1];
      arma::mat un_marg=preobs["un"];
      arma::mat rc_marg=preobs["rc"];
      
      arma::vec obs=un_marg.col(0); // uncensored obs
      arma::vec weight=un_marg.col(1); // uncensored weights
      
      arma::vec rcens=rc_marg.col(0); // right censored obs
      arma::vec rcweight=rc_marg.col(1); // right-censored weights
      
      Zmeanj.col(m)=marginal_expectation(rew, pos, N, alphaj, Sj, obs, weight, rcens, rcweight);
    }
  }
  
  //E-step for E[Z|Ys=ys]
  ZmeanTot=Zmeanj*ed;
  
  // // M step
  for (int i{0}; i < p; ++i) {
    alpha(i) = Bmean(i,0) / (SumOfWeights + SumOfCensored);
    if (alpha(i) < 0) {
      alpha(i) = 0;
    }
    t(i,0) = Nmean(i,p) / Zmean(i,0);
    if (t(i,0) < 0) {
      t(i,0) = 0;
    }
    S(i,i) = -t(i,0);
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
  //
  //additional rj(k)
  
  for (int i=0; i<p; ++i){
    for(int j=0; j<d; ++j){
      
      if(Zmean(i,0)==0 || Zmeanj(i,j)<=0 || ZmeanTot(i,0)==0 ){R(i,j)=0;
      }else{
        R(i,j)=Zmeanj(i,j)/ZmeanTot(i,0);
        
        //R(i,j)=Zmeanj(i,j)/Zmean(i,0);
      }
    }
  }
  rew_sanity_check(R, Rtol); // to avoid numerical instability we consider a reward 0 if it is smaller than Rtol
  
}
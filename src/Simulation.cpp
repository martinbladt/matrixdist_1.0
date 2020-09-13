#include <Rcpp.h>
using namespace Rcpp;

//' Embeded Markov chain of a sub-intensity matrix
//' 
//' Returns the transition probabilities of the embeded Markov chain determined the sub-intensity matrix 
//' @param T A sub-intensity matrix
//' @return The embeded Markov chain
//' @examples
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' embeddedMC(T)
// [[Rcpp::export]]
NumericMatrix embeddedMC(NumericMatrix T, NumericVector t) {
  int p = t.size();
  NumericMatrix Q(p + 1, p + 1);
  
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p + 1; ++j) {
      if (j != i && j < p) {
        Q(i,j) = -1.0 * T(i,j) / T(i,i);
      }
      else if(j == p) {
        Q(i,j) = -1.0 * t[i] / T(i,i);
      }
    }
  }
  Q(p,p) = 1;
  
  return (Q);
}

//' Cumulate matrix
//' 
//' Creates a new matrix with entries the cumulated rows of \code{A}
//' @param A A matrix
//' @return The cumulated matrix
// [[Rcpp::export]]
NumericMatrix cumulateMatrix(NumericMatrix A) {
  int p1 = A.nrow();
  int p2 = A.ncol();
  
  NumericMatrix cumulated(p1, p2);
  
  for (int i = 0; i < p1; ++i) {
    for (int j = 0; j < p2; ++j) {
      if (j == 0){
        cumulated(i,j) = A(i,j);
      }
      else {
        cumulated(i,j) = cumulated(i,j-1) + A(i,j);
      }
    }
  }
  return (cumulated);
}

//' Cumulate vector
//' 
//' Creates a new vector with entries the cumulated entries of \code{A}
//' @param A A vector
//' @return The cumulated vector
// [[Rcpp::export]]
NumericVector cumulateVector(NumericVector A) {
  int p = A.size();
  
  NumericVector cumulated(p);
  
  for (int i = 0; i < p; ++i) {
    if (i == 0){
      cumulated[i] = A[i];
    }
    else {
      cumulated[i] = cumulated[i - 1] + A[i];
    }
  }
  return (cumulated);
}

//' Initial state of Markov jump process
//' 
//' Given the accumulated values of the initial probabilities \code{Pi} and a uniform value \code{u}, it returns the initial state of a Markov jump process
//' @param cumulatedPi A vector
//' @param u A random value in (0,1)
//' @return The initial state of the Markov jump process
// [[Rcpp::export]]
long initialState(NumericVector cumulatedPi, double u) {
  //Given the accumulated values of the initial probabilities (Pi) and a uniform value (u) returns the states that satisfies cumPi_(k-1)<u<cumPi_(k)
  
  if (u <= cumulatedPi[0]) {
    return (0);
  }
  
  for( int i = 1; i < cumulatedPi.size(); ++i) {
    if (cumulatedPi[i - 1] < u && u <= cumulatedPi[i]) {
      return (i);
    }
  }
  return (0);
}

//' New state in a Markov jump process
//' 
//' Given a transition matrix \code{Q}, a uniform value \code{u}, and a previous state \code{k}, it returns the new state of a Markov jump process
//' @param previousState Previous state of the Markov jump process
//' @param cumulatedEmbeddedMC A transition matrix
//' @param u A random value in (0,1)
//' @return The next state of the Markov jump process
// [[Rcpp::export]]
long newState(long previousState, NumericMatrix cumulatedEmbeddedMC, double u) {
  if (u <= cumulatedEmbeddedMC(previousState,0)) {
    return (0);
  }
  
  for (int i = 1; i < cumulatedEmbeddedMC.ncol(); ++i) {
    if (cumulatedEmbeddedMC(previousState,i - 1) < u && u <= cumulatedEmbeddedMC(previousState,i)) {
      return (i);
    }
  }
  
  return (0);
}


//' Random phase-type
//' 
//' Generates a sample of size \code{n} from a phase-type distribution with parameters \code{pi and} \code{T}
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @return The simulated sample
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' t <- -T%*%rep(1, length(T[,1]))
//' n <- 10
//' rphasetype(n, alpha, T, t) 
// [[Rcpp::export]]
NumericVector rphasetype(int n, NumericVector pi, NumericMatrix T, NumericVector t) {
  
  NumericVector sample(n);
  
  NumericMatrix cumulatedEmbeddedMC = cumulateMatrix(embeddedMC(T, t));
  NumericVector cumulatedPi = cumulateVector(pi);
  
  int p = pi.size();
  long state = 0;
  for (int i = 0; i < n; ++i) {
    double time = 0.0;
    state = initialState(cumulatedPi, runif(1)[0]);
    while (state != p) {
      time += log(1.0 - runif(1)[0]) / T(state,state);
      state = newState(state, cumulatedEmbeddedMC, runif(1)[0]);
    }
    sample[i] = time;
  }
  return (sample);
}
#include <Rcpp.h>
using namespace Rcpp;
#include "matrix_functions.h"


//' Embeded Markov chain of a sub-intensity matrix
//' 
//' Returns the transition probabilities of the embeded Markov chain determined the sub-intensity matrix 
//' @param T A sub-intensity matrix
//' @return The embeded Markov chain
//' @examples
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' embeddedMC(T)
// [[Rcpp::export]]
NumericMatrix embeddedMC(NumericMatrix T) {
  long p{T.nrow()};
  NumericMatrix Q(p + 1, p + 1);
  
  NumericVector ee(p, 1);
  NumericMatrix m_e(p, 1, ee.begin());
  NumericMatrix t = matrix_product(T * (-1), m_e);
  
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p + 1; ++j) {
      if (j != i && j < p) {
        Q(i,j) = -1.0 * T(i,j) / T(i,i);
      }
      else if(j == p) {
        Q(i,j) = -1.0 * t(i,0) / T(i,i);
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
//' Generates a sample of size \code{n} from a phase-type distribution with parameters \code{pi} and \code{T}
//' @parm n Sample size
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @return The simulated sample
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' n <- 10
//' rphasetype(n, alpha, T) 
// [[Rcpp::export]]
NumericVector rphasetype(int n, NumericVector pi, NumericMatrix T) {
  
  NumericVector sample(n);
  
  NumericMatrix cumulatedEmbeddedMC = cumulateMatrix(embeddedMC(T));
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



//' Random inhomogeneous phase-type
//' 
//' Generates a sample of size \code{n} from an inhomogeneous phase-type distribution with parameters \code{pi}, \code{T} and \code{beta}
//' @parm n Sample size
//' @parm dist_type Type of IPH: "Pareto", "Weibull", "Gompertz"
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param beta Parameter of the transformation
//' @return The simulated sample
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' beta <- 0.5
//' n <- 10
//' riph(n, "Pareto", alpha, T, beta) 
// [[Rcpp::export]]
NumericVector riph(int n, String dist_type, NumericVector pi, NumericMatrix T, NumericVector beta) {
  
  NumericVector sample(n);
  
  NumericMatrix cumulatedEmbeddedMC = cumulateMatrix(embeddedMC(T));
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
    if (dist_type == "Pareto") {
      time = beta[0] * (exp(time) - 1);
    }
    else if (dist_type == "Weibull") {
      time = pow(time, 1.0 / beta[0]);
    }
    else if (dist_type == "LogLogistic") {
      time = beta[0] * pow(exp(time) - 1, 1 / beta[1]);
    }
    else if (dist_type == "Gompertz") {
      time = log(beta[0] * time + 1) / beta[0];
    }
    sample[i] = time;
  }
  return (sample);
}


//' Random matrix GEVD
//' 
//' Generates a sample of size \code{n} from an inhomogeneous phase-type distribution with parameters \code{pi}, \code{T} and \code{beta}
//' @parm n Sample size
//' @parm dist_type Type of IPH: "Pareto", "Weibull", "Gompertz"
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param mu Location parameter
//' @param sigma Scale parameter
//' @param xi Shape parameter: Default 0 which corresponds to the Gumbel case
//' @return The simulated sample
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' mu <- 3
//' sigma <- 2
//' xi <- 0.5
//' n <- 10
//' rmatrixGEVD(n, alpha, T, mu, sigma, xi) 
//' rmatrixGEVD(n, alpha, T, mu, sigma) 
// [[Rcpp::export]]
NumericVector rmatrixGEVD(int n, NumericVector pi, NumericMatrix T, double mu, double sigma, double xi = 0) {
  
  NumericVector sample(n);
  
  NumericMatrix cumulatedEmbeddedMC = cumulateMatrix(embeddedMC(T));
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
    if (xi == 0) {
      time = mu - sigma * log(time);
    }
    else {
      time = mu + sigma  * (pow(time, - xi) - 1) / xi;
    }
    sample[i] = time;
  }
  return (sample);
}



//' Random inhomogeneous phase-type
//' 
//' Generates a sample of size \code{n} from an inhomogeneous phase-type distribution with parameters \code{pi}, \code{T} and \code{beta}
//' @parm n Sample size
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param beta Parameter of the transformation
//' @return The simulated sample
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' g <- function(x, beta) { x^(1/beta) }
//' beta <- 0.5
//' n <- 10
//' riphfn(n, alpha, T, g, beta) 
// [[Rcpp::export]]
NumericVector riphfn(int n, NumericVector pi, NumericMatrix T, Function g ,NumericVector beta) {
  
  NumericVector sample(n);
  
  NumericMatrix cumulatedEmbeddedMC = cumulateMatrix(embeddedMC(T));
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
    NumericVector g_val = g(time, beta);
    sample[i] = g_val[0];
  }
  return (sample);
}


//' Random MPH*
//' 
//' Generates a sample of size \code{n} from a MPH* distribution with parameters \code{pi}, \code{T} and \code{R}
//' @parm n Sample size
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @return The simulated sample
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' R <- matrix(c(c(1,0,0.8),c(0,1,0.2)), nrow = 3, ncol = 2)
//' n <- 10
//' rmph(n, alpha, T, R) 
// [[Rcpp::export]]
NumericMatrix rmph(int n, NumericVector pi, NumericMatrix T, NumericMatrix R) {
  long dim{R.ncol()};
  
  NumericMatrix sample(n, dim);
  
  NumericMatrix cumulatedEmbeddedMC = cumulateMatrix(embeddedMC(T));
  NumericVector cumulatedPi = cumulateVector(pi);
  
  int p = pi.size();
  long state = 0;
  for (int i = 0; i < n; ++i) {
    double time = 0.0;
    state = initialState(cumulatedPi, runif(1)[0]);
    while (state != p) {
      time = log(1.0 - runif(1)[0]) / T(state,state);
      for (int j{0}; j < dim; ++j) {
        sample(i,j) += R(state, j) * time;
      }
      state = newState(state, cumulatedEmbeddedMC, runif(1)[0]);
    }
  }
  return (sample);
}



//' Random Inhomogeneous MPH*
//' 
//' Generates a sample of size \code{n} from an Inhomogeneous MPH* distribution with parameters \code{pi}, \code{T} and \code{R}
//' @param n Sample size
//' @param dist_type Type of distribution: "Weibull" "Pareto"
//' @param pi Initial probabilities
//' @param T sub-intensity matrix
//' @param beta parameters of the transformations
//' @return The simulated sample
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' T <- matrix(c(c(-1,0,0),c(1,-2,0),c(0,1,-5)), nrow = 3, ncol = 3)
//' R <- matrix(c(c(1,0,0.8),c(0,1,0.2)), nrow = 3, ncol = 2)
//' beta <- c(0.4, 0.7)
//' n <- 10
//' rimph(n, "Weibull", alpha, T, R, beta) 
// [[Rcpp::export]]
NumericMatrix rimph(int n, String dist_type, NumericVector pi, NumericMatrix T, NumericMatrix R, NumericVector beta) {
  long dim{R.ncol()};
  
  NumericMatrix sample(n, dim);
  
  NumericMatrix cumulatedEmbeddedMC = cumulateMatrix(embeddedMC(T));
  NumericVector cumulatedPi = cumulateVector(pi);
  
  int p = pi.size();
  long state = 0;
  for (int i = 0; i < n; ++i) {
    double time = 0.0;
    state = initialState(cumulatedPi, runif(1)[0]);
    while (state != p) {
      time = log(1.0 - runif(1)[0]) / T(state,state);
      for (int j{0}; j < dim; ++j) {
        sample(i,j) += R(state, j) * time;
      }
      state = newState(state, cumulatedEmbeddedMC, runif(1)[0]);
    }
    for (int j{0}; j < dim; ++j) {
      if (dist_type == "Pareto") {
        sample(i,j) = beta[j] * (exp(sample(i,j)) - 1);
      }
      else if (dist_type == "Weibull") {
        sample(i,j) = pow(sample(i,j), 1.0 / beta[j]);
      }
      else if (dist_type == "Gompertz") {
        sample(i,j) = log(beta[j] * sample(i,j) + 1) / beta[j];
      }
    }
  }
  return (sample);
}



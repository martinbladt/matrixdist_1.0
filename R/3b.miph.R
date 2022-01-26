#' Multivariate Inhomogeneous Phase Type distributions
#'
#' Class of objects for multivariate phase-type distributions.
#'
#' @slot name Name of the phase type distribution.
#' @slot gfun A list comprising of the parameters.
#' @slot scale SCale.
#'
#' @return Class object.
#' @export
#'
setClass("miph",
         contains = c("mph"),
         slots = list(
            gfun="list",
            scale="numeric"
         )
)

#' Constructor Function for inhomogeneous multivariate phase-type distributions
#'
#' @param mph An object of class \linkS4class{mph}.
#' @param alpha A probability vector.
#' @param S A list of sub-intensity matrices.
#' @param structure A vector of valid ph structures.
#' @param dimension The dimension of the ph structure (if provided).
#' @param variables Number of marginals.
#' @param gfun Vector of inhomogeneity transform.
#' @param gfun_pars List of parameters for the inhomogeneity functions.
#' @param scale Scale.
#'
#' @return An object of class \linkS4class{iph}.
#' @export
#'
miph<-function(mph=NULL, #object of class mPH
               gfun=NULL, #vector of gfun for each marginal
               gfun_pars=NULL, #List of gfun parameters for each marginal
               alpha=NULL, #vector
               S=NULL, #List
               structure=NULL,
               dimension=3,
               variables=NULL,
               scale=1){
  
  if (all(is.null(c(gfun, gfun_pars)))) {
    stop("input inhomogeneity function and parameters")
  }
  d <- length(gfun)
  
  if(is.null(mph)){
    mph<-mph(alpha = alpha, S = S, structure = structure, dimension = dimension, variables = d)
  }
  
  if(!all(gfun %in% c("pareto", "weibull", "lognormal", "loglogistic", "gompertz", "gev", "identity"))){
    stop("invalid gfun for",which(gfun_check==F))
  }
  
  if (all(gfun %in%  c("pareto", "weibull", "lognormal", "gompertz"))){
    for(i in 1:d){
      if (is.null(gfun_pars[[i]])) gfun_pars[[i]] <- 1
      if (length(gfun_pars[[i]]) != 1 | sum(gfun_pars[[i]] <= 0) > 0) {
        stop(paste("gfun parameter for marginal",i,"should be positive and of length one"))
      } else {
        names(gfun_pars[[i]]) <- paste("beta",i,sep="")
      }
    }
  }
  
  if (any(gfun %in% c("gev"))) {
    index_gev <- which(gfun=="gev")
    for(i in index_gev){
      index_gev2 <- which(gfun[[i]]=="gev")
      if (is.null(gfun_pars[[i]])) gfun_pars[[i]] <- c(0, 1, 1)
      if (length(gfun_pars[[i]]) != 3 | (gfun_pars[[i]][2] > 0) == FALSE) {
        stop("gfun parameter should be of length three: mu, sigma, xi, and sigma > 0")
      } else {
        names(gfun_pars) <- c("mu", "sigma", "xi")
      }
    }
  }
  
  if (any(gfun=="loglogistic")) {
    index_log <- which(gfun=="loglogistic")
    for(i in index_log){
      if (is.null(gfun_pars[[i]])) gfun_pars[[i]] <- c(1, 1)
      if (length(gfun_pars[[i]]) != 2 | (gfun_pars[[i]][1] <= 0) | (gfun_pars[[i]][2] <= 0)) {
        stop("gfun parameter should be positive and of length two: alpha, theta > 0")
      } else {
        names(gfun_pars[[i]]) <- c(paste("alpha",i,sep=""),paste("theta",i,sep=""))
      }
    }
  }
  ginv <-list()
  ginv_prime <- list()
  lambda <- list()
  lambda_prime <- list()
  
  for(i in 1:d){
    f1 <- function(beta, t) t^(beta)
    f2 <- function(beta, t) log(t / beta + 1)
    f3 <- function(beta, t) log(t + 1)^(beta)
    f4 <- function(beta, t) log((t / beta[1])^(beta[2]) + 1)
    f5 <- function(beta, t) (exp(t * beta) - 1) / beta
    f6 <- function(beta, t, w) revers_data_trans(t, w, beta)
    nb <- which(gfun[i] == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
    ginv[[i]] <- base::eval(parse(text = paste("f", nb, sep = "")))
    
    f1 <- function(beta, t) t^(beta) * log(t)
    f2 <- function(beta, t) -t / (beta * t + beta^2)
    f3 <- function(beta, t) log(t + 1)^(beta) * log(log(t + 1))
    f4 <- NA
    f5 <- function(beta, t) exp(t * beta) * (t * beta - 1) / beta^2
    f6 <- NA
    nb <- which(gfun[i] == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
    ginv_prime[[i]] <- base::eval(parse(text = paste("f", nb, sep = "")))
    
    f1 <- function(beta, t) beta * t^(beta - 1)
    f2 <- function(beta, t) (t + beta)^(-1)
    f3 <- function(beta, t) beta * log(t + 1)^(beta - 1) / (t + 1)
    f4 <- NA
    f5 <- function(beta, t) exp(t * beta)
    f6 <- NA
    nb <- which(gfun[i] == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
    lambda[[i]] <- base::eval(parse(text = paste("f", nb, sep = "")))
    
    f1 <- function(beta, t) t^(beta - 1) + beta * t^(beta - 1) * log(t)
    f2 <- function(beta, t) -(t + beta)^(-2)
    f3 <- function(beta, t) log(t + 1)^(beta - 1) / (t + 1) + beta * log(t + 1)^(beta - 1) * log(log(t + 1)) / (t + 1)
    f4 <- NA
    f5 <- function(beta, t) t * exp(t * beta)
    f6 <- NA
    nb <- which(gfun[i] == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
    lambda_prime[[i]] <- base::eval(parse(text = paste("f", nb, sep = "")))
    
  }
  name <- if (is(mph, "miph")) mph@name else paste("inhomogeneous ", mph@name, sep = "")
  
  methods::new("miph",
               name = name,
               pars = mph@pars,
               gfun = list(
                 name = gfun, #a vector
                 pars = gfun_pars, #a list
                 inverse = ginv, #a list
                 inverse_prime = ginv_prime, #a list 
                 intensity = lambda, #a list
                 intensity_prime = lambda_prime #a list
               ),
               scale = scale,
  )
} 

#' Simulation Method for inhomogeneous multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{miph}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed inhomogeneous
#'  multivariate phase-type variables.
#' @export
#'
setMethod("sim", c(x = "miph"), function(x, n = 1000) {
  name <- x@gfun$name
  pars <- x@gfun$pars
  scale <- x@scale
  
  U <- numeric(0)
  for(i in 1:length(name)){
    if (name[i] %in% c("pareto", "weibull", "lognormal", "loglogistic", "gompertz")) {
      U <-cbind(U, scale * riph(n, name[i], x@pars$alpha, x@pars$S[[i]], pars[[i]]))
    }
    if (name[i] %in% c("gev")) {
      U <- cbind(U, scale * rmatrixgev(n, x@pars$alpha, x@pars$S[[i]], pars[[i]][1], pars[[i]][2], pars[[i]][3]))
    }
  }
  return(U)
})
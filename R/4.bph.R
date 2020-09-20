#' Bivariate Phase Type distributions
#'
#' Class of objects for bivariate phase type distributions
#'
#' @slot name name of the phase type distribution
#' @slot pars a list comprising of the parameters.
#'
#' @return
#' @export
#'
setClass("bph",
         contains = c("mph"),
         slots = list(
           blocks = "list"
         )
)

#' Constructor Function for bivariate phase type distributions
#'
#' @param alpha a probability vector.
#' @param S a sub-intensity matrix.
#' @param structure a valid ph structure
#' @param dimensions the dimensions of the ph structures (if provided)
#'
#' @return An object of class \linkS4class{bph}.
#' @export
#'
#' @examples
bph <- function(alpha = NULL, T11 = NULL, T12 = NULL, T22 = NULL, rand_str = FALSE, dimensions = c(3, 3)) {
  if(!rand_str){
    dimensions <- c(dim(T11)[1], dim(T22)[1])
    L <- list(alpha0 = alpha, 
              T11 = T11, 
              T12 = T12, 
              T22 = T22,
              alpha = c(alpha, numeric(length(alpha))),
              S = merge_matrices(T11, T12, T22),
              R = random_phase_BivPH(dimensions[1], dimensions[2])$R
              )
  }else{
    L <- random_phase_BivPH(dimensions[1], dimensions[2])
    names(L) <- c("alpha0", "T11", "T12", "T22", "alpha", "S", "R" )
  }
  new("bph",
      name = paste("bivariate ph(", dimensions[1], ",", dimensions[2], ")", sep = ""),
      pars = list(alpha = L$alpha, S = L$S),
      rew = list(R = L$R),
      blocks = L
  )
}

#' Show Method for bivariate phase type distributions
#'
#' @param x an object of class \linkS4class{bph}.
#' @export
#'
#' @examples
#'
setMethod("show", "bph", function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("number of variables: ", dim(object@rew$R)[2], "\n", sep = "")
  cat("reward matrix: ", "\n", sep = "")
  print(object@rew$R)
})

#' Density Method for bivariate phase type distributions
#'
#' @param x an object of class \linkS4class{bph}.
#' @param y locations
#'
#' @return Density evaluated at locations
#' @export
#'
#' @examples
#'
setMethod("dens", c(x = "bph"), function(x, y = matrix(c(1, 2, 1, 2), 2, 2)) {
  dens <- bivphden(y, x@blocks$alpha0, x@blocks$T11, x@blocks$T12, x@blocks$T22)
  return(list(y = y, dens = dens))
})

#' Distribution Method for bivariate phase type distributions
#'
#' @param x an object of class \linkS4class{bph}.
#' @param y locations
#'
#' @return CDF evaluated at locations
#' @export
#'
#' @examples
#'
setMethod("cdf", c(x = "bph"), function(x, 
                                        q = seq(0, quan(x, .95)$quantile, length.out = 10),
                                        lower.tail = TRUE) {
  #YET TO BE WRITTEN...
  return(list(q = q, cdf = cdf))
})


#' Fit Method for bph Class
#'
#' @param x an object of class \linkS4class{bph}.
#' @param y vector or data.
#'
#' @return An object of class \linkS4class{ph}.
#' @export
#'
#' @examples
#'
setMethod(
  "fit", c(x = "bph", y = "ANY"),
  function(x,
           y,
           weight = numeric(0),
           stepsEM = 1000) {
    is_ibph <- is(x, "ibph")
    if(is_ibph){
      name <- x@gfun$name
      par_g <- x@gfun$pars
      specs <- g_specs(name) 
      inv_g <- specs$inv_g 
      fn <- eval(parse(text = paste("bivm", name, "den", sep = "")))
      mLL <- function(y, alpha, T11, T12, T22, beta) - sum(log(fn(y, alpha, T11, T12, T22, beta)))
    }
    if(length(weight) == 0) weight <- rep(1, dim(y)[1])
     
    ph_par <- x@blocks
    alpha_fit <- clone_vector(ph_par$alpha0)
    T11_fit <- clone_matrix(ph_par$T11)
    T12_fit <- clone_matrix(ph_par$T12)
    T22_fit <- clone_matrix(ph_par$T22)
    
    if(!is_ibph){
      for (k in 1:stepsEM) {
        EMstep_bivph(y, weight, alpha_fit, T11_fit, T12_fit, T22_fit)
        if (k %% 10 == 0) {
          cat("\r", "iteration:", k,
              ", logLik:", sum(log(bivphden(y, alpha_fit, T11_fit, T12_fit, T22_fit))), #weights here??
              sep = " ")
        }
      }
      cat("\n", sep = "")
      z <- bph(alpha_fit, T11_fit, T12_fit, T22_fit)
    }
    if(is_ibph){
      for (k in 1:stepsEM) {
        trans <- cbind(inv_g(y[,1], weight, par_g[1])$obs, inv_g(y[,2], weight, par_g[2])$obs)
        EMstep_bivph(trans, weight, alpha_fit, T11_fit, T12_fit, T22_fit)
        
        opt <- suppressWarnings(optim(par = par_g, fn = mLL,  y = y, alpha = alpha_fit, T11 = T11_fit, T12 = T12_fit, T22 = T22_fit))
        par_g <- opt$par
        if (k %% 10 == 0) {
          cat("\r", "iteration:", k,
              ", logLik:", - opt$value,
              sep = " ")
        }
      }
      cat("\n", sep = "")
      z <- bph(alpha_fit, T11_fit, T12_fit, T22_fit)
      z <- ibph(z, gfun = name, gfun_pars = par_g)
    }
    return(z)
  }
)







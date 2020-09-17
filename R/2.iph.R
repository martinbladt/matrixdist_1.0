#' Inhomogeneous Phase Type distributions
#'
#' Class of objects for inhomogeneous phase type distributions
#'
#' @slot name name of the phase type distribution
#' @slot pars a list comprising of the parameters.
#'
#' @return
#' @export
#'
setClass("iph",
         contains = c("ph"),
         slots = list(
           ph = "ph",
           gfun = "list"
         ),
         prototype = list(
           ph = new("ph"),
           gfun = list()
         )
)

#' Constructor Function for inhomogeneous phase type distributions
#'
#' @param alpha a probability vector.
#' @param S a sub-intensity matrix.
#' @param structure a valid ph structure
#' @param dimension the dimension of the ph structure (if provided)
#' @param fun inhomogeneity transform
#' @param funpars the parameters of the inhomogeneity function
#'
#' @return An object of class \linkS4class{iph}.
#' @export
#'
#' @examples
iph <- function(ph = NULL, gfun = NULL, gfun_pars = NULL, alpha = NULL, S = NULL, structure = NULL, dimension = 3) {
  if(all(is.null(c(gfun, gfun_pars)))){stop("input inhomogeneity function and parameters")}
  if(is.null(A)){ph <- ph(alpha = alpha, S = S, structure = structure, dimension = dimension)}
  if(!gfun %in% c("Pareto", "Weibull", "Gompertz", "GEVD")){stop("invalid gfun")}
  if(gfun %in% c("Pareto", "Weibull", "Gompertz")){
    if(length(gfun_pars) != 1 | sum(gfun_pars <= 0) > 0){
      stop("gfun parameter should be positive and of length one")
    }else{
        names(gfun_pars) <- "beta"
      }
  }
  if(gfun %in% c("GEVD")){
    if(length(gfun_pars) != 3 | (gfun_pars[2] > 0) == FALSE ){
      stop("gfun parameter should be length three: mu, sigma, xi, and sigma > 0")
    }else{
      names(gfun_pars) <- c("mu", "sigma", "xi")
      }
  }
  new("iph",
      ph = ph,
      gfun = list(name = gfun, pars = gfun_pars)
  )
}

#' Show Method for inhomogeneous phase type distributions
#'
#' @param x an object of class \linkS4class{iph}.
#' @export
#'
#' @examples
#'
setMethod("show", "iph", function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  if (length(object@name) > 0) {
    cat("phase-type name: ", object@ph@name, "\n", sep = "")
    cat("parameters: ", "\n", sep = "")
    print(object@ph@pars)
    cat("g-function name: ", object@gfun$name, "\n", sep = "")
    cat("parameters: ","\n", sep = "")
    print(object@gfun$pars)
  } else {
    return()
  }
})

#' Simulation Method for inhomogeneous phase type distributions
#'
#' @param x an object of class \linkS4class{iph}.
#' @param n length of realization.
#'
#' @return A realization of inhomogeneous phase type data
#' @export
#'
#' @examples
#'
setMethod("r", c(x = "iph"), function(x, n = 1000) {
  name <- x@gfun$name
  pars <- x@gfun$pars
  if(name %in% c("Pareto", "Weibull", "Gompertz")){
    U <- riph(n, name, x@ph@pars$alpha, x@ph@pars$S, pars)
  } 
  if(name %in% c("GEVD")){
    U <- rmatrixGEVD(n, x@ph@pars$alpha, x@ph@pars$S, pars[1], pars[2], pars[3])
  }
  return(U)
})

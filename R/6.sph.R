#' Survival Analysis for Phase Type distributions
#'
#' Class of objects for inhomogeneous phase type distributions
#'
#' @slot name name of the phase type distribution
#' @slot pars a list comprising of the parameters.
#'
#' @return
#' @export
#'
setClass("sph",
         contains = c("iph"),
         slots = list(
           coefs = "list",
           type = "character"
         )
)

#' Constructor Function for Survival phase type objects
#'
#' @param ph An object of class \linkS4class{ph}
#' @param iph An object of class \linkS4class{iph}
#'
#' @return An object of class \linkS4class{sph}.
#' @export
#'
#' @examples
sph <- function(ph = NULL, iph = NULL, coefs = list(B = numeric(0), C = numeric(0)), type = "reg") {
  if(!type %in% c("reg", "reg2", "acf")) stop("type must be one of : reg, reg2, acf")
  if (is.null(iph)) {
    if (is.null(ph)) {
      iph <- ph(structure = "General")
    }else{iph <- ph}
    gfun <- list(name = "Indentity", pars = numeric(0))
  }else{
    gfun <- iph@gfun
  }
  iph@pars
  new("sph",
      name = paste("survival", type, iph@name, sep = " "),
      pars = iph@pars,
      gfun = gfun,
      coefs = coefs,
      type = type
  )
}

#' Show Method for survival phase type objects
#'
#' @param x an object of class \linkS4class{sph}.
#' @export
#'
#' @examples
#'
setMethod("show", "sph", function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("g-function name: ", object@gfun$name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@gfun$pars)
  cat("coefficients: ", "\n", sep = "")
  print(object@coefs)
})


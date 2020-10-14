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
           gfun = "list",
           coefs = "list"
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
sph <- function(ph = NULL, iph = NULL, coefs = list(B = NULL, C = NULL)) {
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
      name = paste("survival ", iph@name, sep = ""),
      pars = iph@pars,
      gfun = gfun,
      coefs = coefs
  )
}

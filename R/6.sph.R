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


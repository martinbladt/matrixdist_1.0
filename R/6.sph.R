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
#' @param x An object of class \linkS4class{ph}
#'
#' @return An object of class \linkS4class{sph}.
#' @export
#'
#' @examples
sph <- function(x = NULL, coefs = list(B = numeric(0), C = numeric(0)), type = "reg") {
  if(!type %in% c("reg", "reg2", "aft")) stop("type must be one of : reg, reg2, aft")
  if (!is(x, "iph")) {
    if (!is(x, "ph")) {
      ph <- ph(structure = "General")
    }
    gfun <- list(name = "Identity", pars = numeric(0))
  }else{
    gfun <- x@gfun
  }
  new("sph",
      name = paste("survival", type, x@name, sep = " "),
      pars = x@pars,
      gfun = gfun,
      coefs = coefs,
      type = type,
      scale = 1,
      fit = x@fit
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

#' Coef Method for sph Class
#'
#' @param object an object of class \linkS4class{sph}.
#'
#' @return parameters of sph model
#' @export
#'
setMethod("coef", c(object = "sph"), function(object) {
  L <- append(object@pars, unname(object@gfun$pars))
  names(L)[3:(2 + length(object@gfun$pars))] <- names(object@gfun$pars)
  append(L, object@coefs)
})

#' Evaluation Method for sph Class
#'
#' @param object an object of class \linkS4class{sph}.
#'
#' @return a ph model
#' @export
#'
setMethod("eval", c(x = "sph"), function(x, subject) {
  if(x@gfun$name == "Identity"){
    if(x@type == "aft"){
      z <- ph(alpha = x@pars$alpha,
              S = x@pars$S/as.numeric(exp(subject%*%x@coefs$B))
              )
    }
    else if(x@type == "reg"){
      z <- ph(alpha = x@pars$alpha,
              S = as.numeric(exp(subject%*%x@coefs$B)) * x@pars$S
      )
    }
  }else{
    if(x@type == "aft"){
      z <- iph(gfun = x@gfun$name,
               gfun_pars = x@gfun$pars,
               alpha = x@pars$alpha,
               S = x@pars$S,
               scale = as.numeric(exp(subject%*%x@coefs$B))
      )
    }
    else if(x@type == "reg"){
      z <- iph(gfun = x@gfun$name,
               gfun_pars = x@gfun$pars,
               alpha = x@pars$alpha,
               S = as.numeric(exp(subject%*%x@coefs$B)) * x@pars$S,
      )
    }
    else if(x@type == "reg2"){
      z <- iph(gfun = x@gfun$name,
               gfun_pars = x@gfun$pars * as.numeric(exp(subject%*%x@coefs$C)),
               alpha = x@pars$alpha,
               S = as.numeric(exp(subject%*%x@coefs$B)) * x@pars$S,
      )
    }
  }
  return(z)
})
  

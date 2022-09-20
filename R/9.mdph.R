#' Multivariate Discrete Phase Type distributions
#'
#' Class of objects for multivariate discrete phase-type distributions.
#'
#' @slot name Name of the discrete phase type distribution.
#' @slot pars A list comprising of the parameters.
#' @slot fit A list containing estimation information.
#'
#' @return Class object.
#' @export
#'
setClass("mdph",
         slots = list(
           name = "character",
           pars = "list",
           fit = "list"
         )
)

#' Constructor Function for multivariate discrete phase-type distributions
#'
#' @param alpha A probability vector.
#' @param S A list of sub-transition matrices.
#' @param structure A vector of valid ph structures.
#' @param dimension The dimension of the dph structure (if provided).
#' @param variables The dimension of the multivariate discrete phase-type.
#'
#' @return An object of class \linkS4class{mdph}.
#' @export
#'
mdph <- function(alpha = NULL, S = NULL, structure = NULL, dimension = 3, variables = NULL) {
  if (any(is.null(alpha)) & any(is.null(S)) & is.null(structure)) {
    stop("input a vector and matrix, or a structure")
  }
  if (is.null(variables)) {
    variables <- length(structure)
  }
  if (!any(is.null(structure))) {
    rs <- random_structure(dimension, structure = structure[1])
    alpha <- rs[[1]]
    Sa <- rs[[2]]
    a <- max_diagonal(Sa * (-1)) * (1 + stats::runif(1))
    S <- list()
    S[[1]] <- (a * diag(dimension) + Sa) / a
    for (i in 2:variables) {
      Sa <- random_structure(dimension, structure = structure[i])[[2]]
      a <- max_diagonal(Sa * (-1)) * (1 + stats::runif(1))
      S[[i]] <- (a * diag(dimension) + Sa) / a
    }
    name <- structure
  } else {
    name <- "custom"
  }
  methods::new("mdph",
               name = paste(name, " mdph(", length(alpha), ")", sep = " "),
               pars = list(alpha = alpha, S = S)
  )
}

#' Show Method for multivariate discrete phase-type distributions
#'
#' @param object An object of class \linkS4class{mdph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "mdph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@pars)
  cat("number of variables: ", length(object@pars$S), "\n", sep = "")
})

#' Coef method for mdph class
#'
#' @param object An object of class \linkS4class{mdph}.
#'
#' @return Parameters of multivariate discrete phase-type model.
#' @export
#'
#' @examples
#' obj <- mdph(structure = c("general", "general"))
#' coef(obj)
setMethod("coef", c(object = "mdph"), function(object) {
  object@pars
})

#' Simulation Method for multivariate discrete phase-type distributions
#'
#' @param x An object of class \linkS4class{mdph}.
#' @param n Length of realization.
#' @param equal_marginals Non-negative integer. If positive, it specifies
#' the number of marginals to simulate from, all from the first matrix.
#'
#' @return A realization of a multivariate discrete phase-type distribution.
#' @export
#'
setMethod("sim", c(x = "mdph"), function(x, n = 1000, equal_marginals = 0) {
  
  if(is.vector(x@pars$alpha)) p <- length(x@pars$alpha)
  if(is.matrix(x@pars$alpha)) p <- ncol(x@pars$alpha)
  
  if (equal_marginals == 0) {
    d <- length(x@pars$S)
    
    trans_mat <- list()
    for (j in 1:d) {
      exit_vec <- 1 - rowSums(x@pars$S[[j]])
      t_mat <- cbind(x@pars$S[[j]], exit_vec)
      aux_vec <- rep(0, p + 1)
      aux_vec[p + 1] <- 1
      trans_mat[[j]] <- rbind(t_mat, aux_vec)
    }
    
    states <- 1:p
    result <- matrix(NA, n, d)
    for (i in 1:n) {
      state <- sample(states, 1)
      in_vect <- rep(0, p)
      in_vect[state] <- 1
      for (j in 1:d) {
        result[i, j] <- rdphasetype(1, in_vect, trans_mat[[j]])
      }
    }
  } else {
    d <- equal_marginals
    
    exit_vec <- 1 - rowSums(x@pars$S[[1]])
    trans_mat <- cbind(x@pars$S[[1]], exit_vec)
    aux_vec <- rep(0, p + 1)
    aux_vec[p + 1] <- 1
    trans_mat <- rbind(trans_mat, aux_vec)
    
    states <- 1:p
    result <- matrix(NA, n, d)
    for (i in 1:n) {
      state <- sample(states, 1)
      in_vect <- rep(0, p)
      in_vect[state] <- 1
      for (j in 1:d) {
        result[i, j] <- rdphasetype(1, in_vect, trans_mat)
      }
    }
  }
  return(result)
})

#' Marginal method for mdph class
#'
#' @param x An object of class \linkS4class{mdph}.
#' @param mar Indicator of which marginal.
#' @return An object of the of class \linkS4class{dph}.
#' @export
#'
#' @examples
#' obj <- mdph(structure = c("general", "general"))
#' marginal(obj, 1)
setMethod("marginal", c(x = "mdph"), function(x, mar = 1) {
  x0 <- dph(alpha = x@pars$alpha, S = x@pars$S[[mar]])
  return(x0)
})


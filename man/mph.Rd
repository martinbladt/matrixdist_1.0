% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.mph.R
\name{mph}
\alias{mph}
\title{Constructor function for multivariate phase-type distributions}
\usage{
mph(alpha = NULL, S = NULL, structure = NULL, dimension = 3, variables = NULL)
}
\arguments{
\item{alpha}{A probability vector.}

\item{S}{A list of sub-intensity matrices.}

\item{structure}{A vector of valid ph structures.}

\item{dimension}{The dimension of the ph structure (if provided).}

\item{variables}{The dimension of the multivariate phase-type.}
}
\value{
An object of class \linkS4class{mph}.
}
\description{
Constructor function for multivariate phase-type distributions
}
\examples{
mph(structure = c("gcoxian", "general"), dimension = 5)
}

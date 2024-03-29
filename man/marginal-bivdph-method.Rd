% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8.bivdph.R
\name{marginal,bivdph-method}
\alias{marginal,bivdph-method}
\title{Marginal method for bivdph class}
\usage{
\S4method{marginal}{bivdph}(x, mar = 1)
}
\arguments{
\item{x}{An object of class \linkS4class{bivdph}.}

\item{mar}{Indicator of which marginal.}
}
\value{
An object of the of class \linkS4class{dph}.
}
\description{
Marginal method for bivdph class
}
\examples{
obj <- bivdph(dimensions = c(3, 3))
marginal(obj, 1)
}

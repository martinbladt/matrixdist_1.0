% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/9.mdph.R
\name{marginal,mdph-method}
\alias{marginal,mdph-method}
\title{Marginal method for mdph class}
\usage{
\S4method{marginal}{mdph}(x, mar = 1)
}
\arguments{
\item{x}{An object of class \linkS4class{mdph}.}

\item{mar}{Indicator of which marginal.}
}
\value{
An object of the of class \linkS4class{dph}.
}
\description{
Marginal method for mdph class
}
\examples{
obj <- mdph(structure = c("general", "general"))
marginal(obj, 1)
}

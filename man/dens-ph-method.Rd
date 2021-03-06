% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1.ph.R
\name{dens,ph-method}
\alias{dens,ph-method}
\title{Density Method for phase type distributions}
\usage{
\S4method{dens}{ph}(x, y)
}
\arguments{
\item{x}{an object of class \linkS4class{ph}.}

\item{y}{a vector of locations.}
}
\value{
A list containing the locations and corresponding density evaluations.
}
\description{
Density Method for phase type distributions
}
\examples{
obj <- ph(structure = "general")
dens(obj, c(1, 2, 3))
}

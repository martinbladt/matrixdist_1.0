% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{phdensity}
\alias{phdensity}
\title{Phase-type density}
\usage{
phdensity(x, alpha, S)
}
\arguments{
\item{x}{Non-negative value.}

\item{alpha}{Initial probabilities.}

\item{S}{Sub-intensity matrix.}
}
\value{
The density at \code{x}.
}
\description{
Computes the density of a phase-type distribution with parameters
\code{alpha} and \code{S} at \code{x}.
}

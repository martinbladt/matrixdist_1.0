% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{mloglogisticcdf}
\alias{mloglogisticcdf}
\title{Matrix-loglogistic cdf}
\usage{
mloglogisticcdf(x, alpha, S, beta, lower_tail = TRUE)
}
\arguments{
\item{x}{Non-negative value.}

\item{alpha}{Initial probabilities.}

\item{S}{Sub-intensity matrix.}

\item{beta}{Transformation parameters.}

\item{lower_tail}{Cdf or tail.}
}
\value{
The cdf (tail) at \code{x}.
}
\description{
Computes the cdf (tail) of a matrix-loglogistic distribution with parameters
\code{alpha}, \code{S} and \code{beta} at \code{x}.
}

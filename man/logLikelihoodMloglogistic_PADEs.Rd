% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{logLikelihoodMloglogistic_PADEs}
\alias{logLikelihoodMloglogistic_PADEs}
\title{Loglikelihood of matrix-loglogistic using Pade}
\usage{
logLikelihoodMloglogistic_PADEs(
  h,
  alpha,
  S,
  beta,
  obs,
  weight,
  rcens,
  rcweight,
  scale1,
  scale2
)
}
\arguments{
\item{h}{nuisance parameter}

\item{alpha}{initial probabilities}

\item{S}{sub-intensity}

\item{beta}{in-homogeneity parameter}

\item{obs}{the observations}

\item{weight}{weight of the observations}

\item{rcens}{censored observations}

\item{rcweight}{weight of the censored observations}

\item{scale1}{scale for observations}

\item{scale2}{scale for censored observations}
}
\description{
Loglikelihood for a sample
}

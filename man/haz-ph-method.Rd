% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1.ph.R
\name{haz,ph-method}
\alias{haz,ph-method}
\title{Hazard rate Method for phase type distributions}
\usage{
\S4method{haz}{ph}(x, y = seq(0, quan(x, 0.95)$quantile, length.out = 10))
}
\arguments{
\item{x}{an object of class \linkS4class{ph}.}

\item{y}{locations}
}
\value{
Hazard rate evaluated at locations
}
\description{
Hazard rate Method for phase type distributions
}
\examples{

}
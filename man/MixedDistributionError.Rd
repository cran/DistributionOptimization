% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MixedDistributionError.R
\name{MixedDistributionError}
\alias{MixedDistributionError}
\title{MixedDistributionError}
\usage{
MixedDistributionError(Means, SDs, Weights, Data, rho = 0.5,
  breaks = NULL, Kernels = NULL, ErrorMethod = "chisquare")
}
\arguments{
\item{Means}{Means of the GMM Components}

\item{SDs}{Standard Deviations of the GMM Components}

\item{Weights}{Weights of the GMM Components}

\item{Data}{Empirical Data based on which the GMM is build}

\item{rho}{Ratio of OverlappingError vs Fitting Error}

\item{breaks}{vector containing the breaks between bins}

\item{Kernels}{positions at which density is to be compared}

\item{ErrorMethod}{"pdeerror": fitting error is measured through Pareto Density Estimation.
"chisquare": fitting error is measured through the Chi Square fitting error.}
}
\value{
Mixed Error
}
\description{
Calculates a fitting error as well as the overlapping measure for the mixtures. Combines them with ratio rho in favor of Overlapping.
}
\examples{
Data = c(rnorm(50,1,2), rnorm(50,3,4))
MixedDistributionError(c(1,3), c(2,4), c(0.5,0.5), Data = Data)
}
\author{
Florian Lerch
}

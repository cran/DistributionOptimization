% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OverlapErrorByDensity.R
\name{OverlapErrorByDensity}
\alias{OverlapErrorByDensity}
\title{OverlapErrorByDensity}
\usage{
OverlapErrorByDensity(Means, SDs, Weights, Data = NULL, Kernels = NULL)
}
\arguments{
\item{Means}{Means of the GMM Components}

\item{SDs}{Standard Deviations of the GMM Components}

\item{Weights}{Weights of the GMM Components}

\item{Data}{Dataset that the GMM should be compared with}

\item{Kernels}{if length(Kernels) = 1: amount of kernels
if length(Kernels) > 1: kernels in dataspace at which the GMM Components will be compared with each other}
}
\value{
List:
OverlapError    Error for estimating the maximal Overlap of AUC of PDFs of each pair of GMM Components
Kernels         Kernels that were used for comparing the GMM Components
}
\description{
Similarity in GMM by Density
}
\details{
Calculates the similarity (overlap) between multiple modes in Gaussian Mixture Models. Kernels at equally distanced positions
are used, if not explicitly given.
}
\examples{
Data = c(rnorm(50,1,2), rnorm(50,3,4))
V<-OverlapErrorByDensity(c(1,3), c(2,4), c(0.5,0.5), Data = Data)
AdaptGauss::PlotMixtures(Data, c(1,3), c(2,4), SingleGausses = TRUE)
print(V$OverlapError)
}
\author{
Florian Lerch
}

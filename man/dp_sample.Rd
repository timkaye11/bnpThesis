% Generated by roxygen2 (4.0.2): do not edit by hand
\name{dp_sample}
\alias{dp_sample}
\title{Sample for the Dirichlet Process}
\usage{
dp_sample(n, alpha, base)
}
\arguments{
\item{n}{The number of samples to sample}

\item{alpha}{The concentration parameter for the dirichlet process}

\item{base}{The base distribution (commonly refered to as H) for the dirichlet process}
}
\description{
Returns a sample (non-unique) from the dirichlet process with a specified
concentration parameter (alpha) and a base distribution
}
\examples{
base <- function(n) rpois(n, 20)
samp <- dp_sample(10, 5, base)
}
\keyword{dirichlet}
\keyword{process}


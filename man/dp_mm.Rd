% Generated by roxygen2 (4.0.2): do not edit by hand
\name{dp_mm}
\alias{dp_mm}
\title{Dirichlet Process Mixture Modeling with Gibbs Sampling}
\usage{
dp_mm(data, alpha_0, num_iters = 100, burnin = 50, useSM = FALSE, ...)
}
\arguments{
\item{data}{The data (a data.frame) to be passed in to the DP mixture model}

\item{alpha_0}{The alpha value (concentration parameter) for the DP}

\item{num_iters}{the number of iterations to run the gibbs sampler}

\item{burnin}{The number of iterations to 'burn' for the Markov Chain}

\item{useSM}{Use the softmax trick to configure the 'best' alpha value. Default is FALSE}
}
\description{
Implements a Dirichlet Process Mixture Model with Gibbs Sampling. The base distribution for the Dirichlet
Process is a Normal-Gamma distribution, which results in a predictive posterior of the t-distribution. After
Gibbs Sampling, this function returns a cluster parameters for each of the clusters, as the Dirichlet Process Mixture Model
is a non-parametric clustering algorithm. The alpha (concentration parameter for the DP) is resampled after each iteration
using the soft-max trick, and each of the cluster parameters are resampled at each iteration as well.
}
\examples{
# Makeup some normal mixture model
x <- c(rnorm(50, 3, 1), rnorm(50, 6, .5), rnorm(50, 9, 1.5))
# hist(x, 10)
output <- dp_mm(x, alpha_0=5,num_iters=100)
plot(output)
summary(output)
}
\keyword{clustering}
\keyword{dirichlet}
\keyword{mixture}
\keyword{model,}
\keyword{process,}


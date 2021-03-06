% Generated by roxygen2 (4.0.2): do not edit by hand
\name{mondrian_1d}
\alias{mondrian_1d}
\title{A R function for a one-dimensional mondrian process}
\usage{
mondrian_1d(budget, interval = c(0, 1))
}
\arguments{
\item{budget}{The budget value for the mondrian process. Higher values lead to more splits}

\item{interval}{The interval to create the splits on. The default is set to [0,1]}
}
\description{
The mondrian process (http://danroy.org/papers/RoyTeh-NIPS-2009.pdf) is a point process
that is used in the Mondrian forest algorithm. It creates partitions based on a budget value,
and a specified interval
}
\examples{
mondrian <- mondrian_1d(5)
mondrian
}
\keyword{mondrian}
\keyword{process}
\keyword{process,}


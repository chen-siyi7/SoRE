\name{pl_logmass}
\alias{pl_logmass}
\title{Plackett-Luce log mass function}
\usage{pl_logmass(r, omega)}
\arguments{
\item{r}{Integer permutation vector of length p.}
\item{omega}{Numeric vector of strictly positive weights, length p.}
}
\value{Scalar log-density.}
\description{Evaluates log PL(r | omega) for a given permutation.}
\examples{
omega <- c(10, 5, 3, 1, 1)
pl_logmass(c(1, 2, 3, 4, 5), omega)
pl_logmass(c(5, 4, 3, 2, 1), omega)
}

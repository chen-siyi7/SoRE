\name{pl_sample}
\alias{pl_sample}
\title{Sample a permutation from the Plackett-Luce distribution}
\usage{pl_sample(omega, n = 1)}
\arguments{
\item{omega}{Numeric vector of strictly positive Plackett-Luce weights.}
\item{n}{Number of permutations to draw (default 1).}
}
\value{Integer vector (n=1) or n-by-p matrix (n>1) of permutations.}
\description{
Draws permutations from PL(omega) using the Gumbel-max construction.
Cost is O(p log p) per draw.
}
\examples{
omega <- c(10, 5, 3, 1, 1)
pl_sample(omega)
table(pl_sample(omega, n = 1000)[, 1]) / 1000
}

\name{build_anchor_weights}
\alias{build_anchor_weights}
\title{Build SoRE anchor weights from functional annotation sources}
\usage{
build_anchor_weights(coding = NULL, chrom_open = NULL,
                     conserved = NULL, z_prior = NULL, p = NULL)
}
\arguments{
\item{coding}{Optional integer vector of coding category (3=LoF, 2=missense, 1=other).}
\item{chrom_open}{Optional logical vector indicating chromatin peak overlap.}
\item{conserved}{Optional logical vector indicating conservation (e.g., GERP > 4).}
\item{z_prior}{Optional numeric vector of prior-GWAS z-scores.}
\item{p}{Optional explicit length when all inputs are NULL.}
}
\value{Numeric vector of strictly positive anchor weights.}
\description{
Constructs the multiplicative anchor weight vector tilde-omega_j used
as input to V-SoRE or SoRE Gibbs.
}
\examples{
p <- 100
w <- build_anchor_weights(
  coding     = sample(1:3, p, replace = TRUE),
  chrom_open = runif(p) < 0.05,
  conserved  = runif(p) < 0.10,
  z_prior    = rnorm(p, sd = 0.5)
)
summary(w)
}

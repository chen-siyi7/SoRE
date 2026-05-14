\name{sore_gibbs}
\alias{sore_gibbs}
\title{Blocked Gibbs sampler for the full SoRE posterior}
\usage{
sore_gibbs(z, R, anchor_weights, n_iter = 10000, burn_in = 5000,
           n_chains = 4, W = NULL, a_tau = 1, b_tau = 1, seed = NULL)
}
\arguments{
\item{z}{Numeric vector of marginal z-scores, length p.}
\item{R}{p-by-p ridge-regularized LD matrix.}
\item{anchor_weights}{Numeric vector of strictly positive anchor weights.}
\item{n_iter}{Number of Gibbs sweeps after burn-in (default 10000).}
\item{burn_in}{Burn-in length (default 5000).}
\item{n_chains}{Number of independent chains (default 4).}
\item{W}{Slab variance (scalar).}
\item{a_tau, b_tau}{Inverse-gamma prior hyperparameters.}
\item{seed}{Optional integer seed.}
}
\value{List with posterior summaries.}
\description{
Blocked Gibbs sampler for the full SoRE posterior. Implements the six
full conditionals described in Web Appendix S4.
}
\seealso{\code{\link{vsore}}}

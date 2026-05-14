\name{credible_sets}
\alias{credible_sets}
\title{Construct credible sets from posterior inclusion probabilities}
\usage{
credible_sets(pip, coverage = 0.95, min_purity = NULL, R = NULL)
}
\arguments{
\item{pip}{Numeric vector of posterior inclusion probabilities, length p.}
\item{coverage}{Nominal credible-set coverage (default 0.95).}
\item{min_purity}{Optional minimum purity (median minimum |LD| within CS).}
\item{R}{Optional p-by-p LD matrix (required if min_purity is set).}
}
\value{List with elements \code{cs}, \code{size}, \code{coverage}, \code{purity}, \code{lead_pip}.}
\description{
Builds credible sets greedily by sorting variants in decreasing PIP order
and adding them until the cumulative sum first exceeds the nominal level.
}
\examples{
pips <- c(0.6, 0.3, 0.05, 0.04, 0.01)
credible_sets(pips, coverage = 0.95)
}

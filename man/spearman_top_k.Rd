\name{spearman_top_k}
\alias{spearman_top_k}
\title{Spearman rank correlation on the top-K subset}
\usage{spearman_top_k(x, y)}
\arguments{\item{x}{Numeric vector of weights or scores.}\item{y}{Numeric vector of PIPs.}}
\value{Scalar Spearman correlation in [-1, 1].}
\description{
Computes Spearman's rank correlation between two vectors restricted to
their top-K positions by the first vector.
}

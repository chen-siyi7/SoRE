\name{tau2_diagnostic}
\alias{tau2_diagnostic}
\title{Interpret the V-SoRE tau-squared diagnostic}
\usage{tau2_diagnostic(tau2)}
\arguments{\item{tau2}{Numeric scalar or vector of tau-squared values.}}
\value{Character vector with values "strong", "moderate", "uninformative", or "no annotation".}
\description{
Translates the estimated tau-squared from a V-SoRE fit into a categorical
interpretation following the thresholds calibrated in the manuscript.
}
\examples{
tau2_diagnostic(c(0.1, 0.5, 1.3, NA))
}

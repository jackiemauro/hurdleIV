% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglik_lgiv.R
\name{loglik_lgiv}
\alias{loglik_lgiv}
\title{Lognormal IV likelihood}
\usage{
loglik_lgiv(t)
}
\arguments{
\item{t}{a vector of parameters to be maximized. Order should be:
covariance matrix elements (leaving out leading 1), then betas, then
gammas, then pis.}
}
\value{
returns the loglikelihood of the given parameter vector for the
dataset.
}
\description{
Returns the likelihood for a lognormal IV. Requires certain
parameters to be attached/in the environment. Can't figure out how to
just pass them to the function without the function later trying to
optimize for them. This should essentially always be used in a wrapper
for that reason.
}

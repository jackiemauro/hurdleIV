% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglik_cragg.R
\name{loglik_craggiv}
\alias{loglik_craggiv}
\title{Cragg IV likelihood}
\usage{
loglik_craggiv(t)
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
Returns the likelihood for a Cragg IV. Requires certain
parameters to be attached/in the environment. Note that this uses the
MG method for the errors, described in cragg_errs. This is following
the method worked out with Max.
}

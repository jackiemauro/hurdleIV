% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapperFunc.R
\name{hurdle.IV}
\alias{hurdle.IV}
\title{Hurdle IV regression}
\usage{
hurdle.IV(formula, inst, endog, exog, data, endog_reg = list(),
  start_val = list(), type = "lognormal", options = list(cholesky = T,
  maxit = 5000, trace = 0, method = "BFGS"))
}
\arguments{
\item{formula}{the second stage regression: y~exogenous + endogenous}

\item{inst}{a vector of your instrument(s): c(inst1,inst2)}

\item{endog}{a vector of your endogenous variable(s): c(end1,end2)}

\item{exog}{a vector of your exogenous variable(s): c(ex1,ex2)}

\item{data}{the dataframe}

\item{endog_reg}{a list of endogenous regression formulae. By default will
estimate endogenous ~ all exogenous variables + all instruments}

\item{start_val}{an optional list of start values. By default the function
will find start values using simple linear/probit regressions of the specified
formulae}

\item{type}{either "lognormal" or "cragg"}

\item{options:}{cholesky true or false; maxit = maximal number of iterations;
trace = 0; method = "BFGS". Options are similar to those for optim, see
optim documentation for more details}
}
\value{
Returns estimated parameters as well as a the hessian and standard
deviations
}
\description{
Run a hurdle IV regression, either lognormal IV or cragg IV.
You should specify the model you want to fit. Currently this does not have
a nice summary function for the output.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_df.R
\name{make.df}
\alias{make.df}
\title{Make a dataset of Normals}
\usage{
make.df(mean, sd, n, pref = F)
}
\arguments{
\item{mean}{a vector of means for each column of your dataset}

\item{sd}{a vector of standard deviations for each column of your dataset}

\item{n}{the length of your variables}

\item{pref}{default is False but fill in if you want to specify the names
of your columns}
}
\value{
returns a dataset of normally distributed variables
}
\description{
Generates a dataset with normally distributed variables.
You can specify the mean and variance of each column.
}

#' Resample for truncated normal errors
#'
#' @description The truncated multivariate normal errors can be generated
#' in one of two ways. Either we sample the first stage and probit errors
#' first and then resample the linear regression's errors until we have
#' positive values, or we resample the vector of all three errors until
#' the linear value is positive. This is the second method, use cragg_errs1
#' for the second, which is less intuitive.
#'
#'
#' @param cov the covariance matrix. This should be untransformed, the
#' terms will be multiplied by the coefficients within the resampling
#' procedure.
#' @param x1 your exogenous variables (a dataframe)
#' @param z your instrument (a dataframe)
#' @param pi a vector of coefficients for the first stage regression
#' @param gamma a vector of coefficients for the second stage probit
#' @param beta a vector of coefficients for the second stage linear regression
#' @param n the number of errors to be generated
#'
#' @return returns a list of your errors and the three generated variables:
#' the endogenous regressor, the censoring variable and the outcome variable
#'
cragg_errs2<-function(cov,pi,x1,gamma,beta,n,z){
  require("MASS")
  j=1
  endog = c(rep(NA,n))
  y0 = c(rep(NA,n))
  yStar = c(rep(NA,n))
  errors = matrix(c(rep(0,n*dim(cov)[1])),ncol = dim(cov)[1])

  while(j<=n){
    err = mvrnorm(1,rep(0,dim(cov)[1]),cov)
    frame = as.matrix(cbind(1,x1[j,],z[j,]))
    endog[j] = frame%*%pi + err[3]
    frame2 = as.matrix(cbind(1,x1[j,],endog[j]))
    y0[j] = as.numeric(frame2%*%gamma + err[1] > 0)
    yStar[j] = frame2%*%beta + err[2]

    if(yStar[j]>0){
      errors[j,] = err
      j = j+1
    }

  }

  return(list(errors = errors, endog = endog, y0 = y0, yStar = yStar))
}




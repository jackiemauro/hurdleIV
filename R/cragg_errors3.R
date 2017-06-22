#' Resample for truncated normal errors
#'
#' @description The truncated multivariate normal errors can be generated
#' in one of three ways. 1) we sample the first stage and probit errors
#' first and then resample the linear regression's errors until we have
#' positive values. 2) we resample the vector of all three errors until
#' the linear value is positive. 3) we send all our variables to zero first
#' then resample the full vector, accepting if it lies in Q1.
#' This is the third method.
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
cragg_errs3<-function(cov,pi,x1,gamma,beta,n,z){
  require("MASS")

  # generate initial dataset
  errors = mvrnorm(n,rep(0,dim(cov)[1]),cov)
  frame = as.matrix(cbind(1,x1,z))
  endog = frame%*%pi + errors[,3]
  frame = as.matrix(cbind(1,x1,endog))
  y0Star =frame%*%gamma + errors[,1]
  y1Star = frame%*%beta + errors[,2]

  # send to zero if y0Star<=0
  y0cut = which(y0Star<=0)
  x2plus = endog; x2plus[y0cut] = 0
  y1plus = y1Star; y1plus[y0cut] = 0

  # resample observations with y0Star > 0 and y1Star < 0
  # only accept if new observation has y0Star > 0 and y1Star > 0
  for(j in which(y1plus<0)){
    while(y1plus[j]<=0){
      err = mvrnorm(1,rep(0,dim(cov)[1]),cov)

      frame = as.matrix(cbind(1,x1[j,],z[j,]))
      a = frame%*%pi + err[3]

      frame2 = as.matrix(cbind(1,x1[j,],endog[j]))
      b = frame2%*%gamma + err[1]
      c = frame2%*%beta + err[2]

      if((c>0) & (b>0)){
        x2plus[j] = a;  y0Star[j] = b; y1plus[j] = c; errors[j,]=err
      }
    }
  }
  return(list(errors = errors, endog = x2plus, y0 = as.numeric(y0Star>0), yStar = y1plus))
}



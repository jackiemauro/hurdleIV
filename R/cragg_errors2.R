#' Resample only linear errors for truncated normal errors
#'
#' @description The truncated multivariate normal errors can be generated
#' in one of two ways. Either we sample the first stage and probit errors
#' first and then resample the linear regression's errors until we have
#' positive values, or we resample the vector of all three errors until
#' the linear value is positive. This is the first method, use cragg_errors
#' for the second, which is more intuitive.
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
  newcov = matrix(c(cov[1,1],cov[1,3],cov[1,3],cov[3,3]), ncol = 2, byrow = T)
  err1 = mvrnorm(n,rep(0,2),newcov)
  frame = as.matrix(cbind(1,x1,z))
  endog = frame%*%pi + err1[,2]
  frame2 = as.matrix(cbind(1,x1,endog))
  y0 = as.numeric(frame2%*%gamma + err1[,1] > 0)

  # get conditional mean and variance
  mu = c(cov[1,2], cov[2,3])%*%solve(newcov)%*%t(err1)
  sig2 = cov[2,2] - c(cov[1,2], cov[2,3])%*%solve(newcov)%*%c(cov[2,3], cov[1,2])

  j=1
  yStar = c(rep(NA,n))
  err2 = c(rep(NA,n))

  while(j<=n){
    err = rnorm(1,mu,sqrt(sig2))
    yStar[j] = frame2[j,]%*%beta + err
    if(yStar[j]>0){
      err2[j] = err
      j = j+1
    }

  }
  errors = cbind(err1,err2)
  return(list(errors = errors, endog = endog, y0 = y0, yStar = yStar))
}

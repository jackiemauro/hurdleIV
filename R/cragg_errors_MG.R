#' Resample for truncated normal errors
#'
#' @description There are a bunch of ways of thinking about this DGP. This is
#' what Max G'Sell had in mind. We generate x2 and y0 first, and then sample
#' from y1 conditionally on those two and on itself being positive.
#' This can only handle one endogenous variable for now. Testing it with
#' more than one exogenous variable.
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
cragg_errs_MG<-function(cov,pi,x1,gamma,beta,n,z){
  require("MASS")
  errs = mvrnorm(n,c(0,0,0),cov)

  #Transformation matrix from (eta, u, v) to (y0*, y1*, x2)
  A = rbind(
    c(1,0,tail(gamma,1)),
    c(0,1,tail(beta,1)),
    c(0,0,1)
  )

  #Covariance of (y0*, y1*, x2)
  Sig = A%*%cov%*%t(A)
  Sig_02 = Sig[-2,-2]

  frame1 = as.matrix(cbind(1,x1,z))
  mu_x2 = frame1%*%pi
  frame2 = as.matrix(cbind(1,x1,mu_x2))
  mu_y1 = frame2%*%beta
  mu_y0 = frame2%*%gamma

  x2 = mu_x2 + errs[,3]
  frame3 = as.matrix(cbind(1,x1,x2))
  y0 = frame3%*%gamma + errs[,1]

  mu_y1_y0x2 = mu_y1 + t(Sig[2,c(1,3),drop=FALSE]%*%solve(Sig_02)%*%t(cbind(y0-mu_y0,x2-mu_x2)))
  sig2_y1_y0x2 = Sig[2,2] - Sig[2,c(1,3),drop=FALSE]%*%solve(Sig_02)%*%Sig[c(1,3),2,drop=FALSE]

  y1 = rep(0,n)
  for(j in which(y0>0)){
    while(y1[j]==0){
      y=rnorm(1,mu_y1_y0x2,sd=sqrt(sig2_y1_y0x2))
      if( y>0 ){y1[j]=y}
    }
  }

  return(list(errors = errs, endog = x2, y0 = as.numeric(y0>0), yStar = y1))
}


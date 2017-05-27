#' Resample for truncated normal errors
#'
#' @description There are a bunch of ways of thinking about this DGP. This is
#' (I think) what Dave Choi had in mind. This can only handle one endogenous
#' variable and one exogenous variable for now. Note that this is very
#' similar to cragg3, but we resample all of y1 when y0 is not zero. Also note
#' that this is a little weird in that the x2 we end up observing is not the
#' exact same one we used in generating y0. So there may be rows in the data
#' where the observed x2 would lead to a negative y0 but we still observe a
#' positive value of y1.
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
cragg_errs_DC<-function(cov,pi,x1,gamma,beta,n,z){
  require("MASS")
  errs = mvrnorm(n,c(0,0,0),cov)

  x2 = pi[1] + pi[2]*x1 + pi[3]*z + errs[,3]
  y0 = gamma[1] +gamma[2]*x1 + gamma[3]*x2 + errs[,1]

  yr3 = as.numeric(y0>=0)

  # redraw all three (call the vector e') as long as y0<0 from P(e'|y0'>0,y1'>0)
  for(j in which(yr3!=0)){
    while(yr3[j]<0){
      err = mvrnorm(1,c(0,0,0),Sig_err)
      t2 = pi[1] + pi[2]*x1 + pi[3]*z + err[3]
      t0 = gamma[1] + gamma[2]*x1 + gamma[3]*t2 + err[1]
      t1 = beta[1] + beta[2]*x1 + beta[3]*t2 + err[2]
      if( (t1>0) & (t0>0) ){yr3[j]=t1; x2[j] = t2;errs[j,]=err}
    }
  }
  return(list(errors = errs, endog = x2, yStar = yr3))
}


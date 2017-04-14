#' Simulate a dataset
#'
#' @description Simulates a dataset from either Cragg or Lognormal distributions.
#' Allows you to specify the data generating process, though it still
#' can't handle more than one endogenous regression
#'
#' @param formula For now, must be False. In the future, will allow you to leave
#' out exogenous variables you use in first stage regressions from second
#' stage regressions.
#' @param pi a vector of the coefficients on your first stage regression.
#' Should be in the order: intercept, exogenous variables, instruments.
#' Defaults to (1,-1,3).
#' @param gamma a vector of the coefficients on your second stage probit
#' regression. Should be in the order: intercept, exogenous variables,
#' endogenous variable. Defaults to (-.2,.8,.07)
#' @param beta a vector of the coefficients on your second stage linear
#' regression. Should be in the order: intercept, exogenous variables,
#' endogenous variable. Defaults to (.05,.06,.02)
#' @param endog_reg for now, leave as an empty list. Will create one endogenous
#' regression which includes the instrument and all exogenous variables
#' @param exog_mean the mean(s) of your exogenous variable(s). Defaults to 1
#' @param exog_sd the standard deviation(s) of your exogenous variable(s).
#' Defaults to 1.
#' @param z_mean the mean(s) of your instrument(s). Defaults to 3.
#' @param z_sd the standard deviation(s) of your instrument(s).
#' Defaults to 1.
#' @param y_sd the standard deviation of your second stage linear regression.
#' Defaults to 2.
#' @param rho the covariance term between second stage linear and probit
#' regressions. Defaults to 0.2.
#' @param tau0 the covariance term between the second stage probit and the
#' first stage regression. Defaults to 0.3.
#' @param tau1 the covariance term between the second stage linear regression and the
#' first stage regression. Defaults to 0.1.
#' @param n the number of observations. Defaults to 10,000.
#' @param type defaults to lognormal. Enter "cragg" for cragg model.
#' @param options Options are: silent = F (will print histograms and percent censored),
#' cragg_errors = {1,2} cragg error type 1 resamples entire vector of errors.
#' See documentation for cragg_errs for more detail.
#'
#' @return returns a simulated dataset

hurdle.IV.sim <- function(formula = F,
                                   pi = c(1,-1,3), #intercept, exog, inst
                                   gamma = c(-.2,.8,.07), #intercept, exog, endog
                                   beta = c(.05,.06,.02), #intercept, exog, endog
                                   endog_reg = list(),
                                   exog_mean = 1,
                                   exog_sd = 1,
                                   z_mean = 3,
                                   z_sd = 1,
                                   endog_sd = 5,
                                   y_sd = 2,
                                   rho = .2,
                                   tau0 = .3,
                                   tau1 = .1,
                                   n = 10000,
                                   options = list(silent = F, cragg_errors = 1),
                                   type = "lognormal"){

  #checks
  n_z = length(z_mean)
  n_x2 = length(endog_sd)

  if(n_x2 > n_z){
    stop("Error: more endogenous variables than instruments")
  }
  if(length(gamma)!= length(beta) ){
    stop("Error: coefficient vector lengths differ")
  }
  if(length(tau0)!=n_x2|length(tau1)!=n_x2){
    stop("Error: length of tau vectors must equal number of endogenous variables")
  }
  if(length(endog_reg)!=0){
    # need to allow you to exclude some exogenous variables
    stop("This functionality not developed yet, set endog_reg = list()")
  }
  if(formula!=F){
    # need to allow you to exclude some exogenous variables
    stop("This functionality not developed yet, set formula = F")
  }

  # make instruments and exogenous variables as random normals
  z = make.df(mean = z_mean, sd = z_sd, n = n, pref = "inst")
  x1 = make.df(mean = exog_mean, sd = exog_sd, n = n, pref = "exog")
  df = data.frame(x1,z)

  # make covariance matrix
  require(MASS)
  cov = make.cov(rho=rho
                 ,tau0=tau0
                 ,tau1=tau1
                 ,y_sd=y_sd
                 ,endog_sd=endog_sd)

  if(type == "cragg"){
    if(options$cragg_errors == 1){
      temp = cragg_errs2(cov=cov,df=df,pi=pi,x1=x1,gamma=gamma,beta=beta,n=n,z=z)
    }
    else{
      temp = cragg_errs(cov=cov,df=df,pi=pi,x1=x1,gamma=gamma,beta=beta,n=n,z=z)
    }
    endog = temp['endog'][[1]]
    errors = temp['errors'][[1]]
    y0 = temp['y0'][[1]]
    yStar = temp['yStar'][[1]]
    y = y0*yStar
  }
  else{
    errors = mvrnorm(n,rep(0,dim(cov)[1]),cov)
    frame = as.matrix(cbind(1,df))
    endog = frame%*%pi + errors[,3]
    frame = as.matrix(cbind(1,x1,endog))
    y0 = as.numeric(frame%*%gamma + errors[,1] > 0)
    yStar = frame%*%beta + errors[,2]
    y = exp(yStar)*y0
  }

  out = data.frame(y,y0,endog,df)

  if(silent == F){
    par(mfrow = c(1,2))
    hist(log(out$y), main = "log(y)",xlab = "")
    hist(out$y, main = "y", xlab = "")
    par(mfrow = c(1,1))
    print(paste(mean(out$y0)*100," percent uncensored"))
  }
  return(out)
}


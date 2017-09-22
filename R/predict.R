#' Predict from lognormal IV
#'
#' @description After getting the output of a lognormal IV, get predicted values
#' for log(y) and y0Star.
#'
#' @param res the output of a lognormal IV regression
#' @param x the covariates (not including intercept)  you want to use in prediction
#'
#' @example
#' dat = hurdle.IV.sim()
#' out = hurdle.IV(y~exog1 + endog, inst = inst1, endog = endog, exog = exog1, data = dat)
#' pred = hurdle.IV.predict(out, dat[,3:4])
#'
#' @return a dataframe with predicted log(y) and y0


hurdle.IV.predict<- function(res, x){
  # input the result of a hurdleIV (lognormal) regression and get a vector of
  # predicted ln(y)'s and y0Star's

  n <- dim(x)[1]
  beta <- res$parameters$beta
  gamma <- res$parameters$gamma

  par.names <- names(beta)[-1]
  cov.names <- names(x)
  if(any(par.names != cov.names)){
    cat("Warning: parameter names different from covariate names",
                "\nParameter names:", par.names, "\nCovariate names:",
                names(x), sep = " ")}

  lny <- c(t(as.matrix(beta)) %*% t(as.matrix(cbind(rep(1,n),x))))
  y0 <- c(t(as.matrix(gamma)) %*% t(as.matrix(cbind(rep(1,n),x))))
  y <- exp(lny)*as.numeric(y0>0)

  return(data.frame(lny.pred = lny, y0.pred = y0, y.pred = y))
}

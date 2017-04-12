#' Get the starting points for your regression
#' 
#' @description The maximum likelihood procedure needs a place to start.
#' This runs traditional IV probit/linear regressions to get starting points.
#'
#' @param formula the second stage linear regression (outcome ~ exogenous + endogenous)
#' @param endog_reg a list of endogenous regressions. By default it will 
#' include all exogenous regressors and one instrument
#' @param data your data, as a dataset
#' @param type either "lognormal" or "cragg"
#' 
#' @return returns a list of starting parameters for the coefficients
#' and the covariance matrix

start.val <- function(formula
                      ,endog_reg
                      ,data
                      ,type){
  
  # this would be better if it used the predicted values of endogenous vars
  # requires updating formulae in a way that escapes me rn
  
  endogs = lapply(endog_reg, function(x) lm(x, data = data))
  endog.preds = lapply(endogs, function(x) predict(x))
  pis = lapply(endogs, function(x) coef(x))
  endog_sd = lapply(endogs, function(x) sd(x$residuals))
  
  if(type == "lognormal"){
    lin.form = update(formula, log(.)~.)
  }
  else if(type == "cragg"){
    lin.form = formula
  }
  
  Iy = data$outcome>0
  data$y0 = as.numeric(Iy)
  linear = lm(lin.form, data = data[Iy,])
  prob.form = update(formula, y0 ~ . ) 
  probit = glm(prob.form, family = binomial("probit"), data = data)
  betas = coef(linear)
  y_sd = sd(linear$residuals)
  gammas = coef(probit)
  
  rho = cov(linear$residuals, probit$residuals[Iy])
  tau0 = lapply(endogs, function(x) cov(probit$residuals, x$residuals))
  tau1 = lapply(endogs, function(x) cov(linear$residuals, x$residuals[Iy]))
  
  return(list(beta = betas,
              gamma = gammas,
              pi = pis,
              endog_sd = unlist(endog_sd),
              y_sd = y_sd,
              tau0 = tau0,
              tau1 = tau1,
              rho = rho))
}
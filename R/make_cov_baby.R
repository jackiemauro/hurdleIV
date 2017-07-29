#' Helper function to make covariance matrix
#'
#' @description Specify elements of the covariance matrix and it organizes
#' them for you
#'
#' @param rho covariance between second stage regressions
#' @param tau0 covariance between second stage probit and first stage
#' regression
#' @param tau1 covariance between second stage linear and first stage
#' regression
#' @param y_sd second stage linear regression standard deviation
#' @param endog_sd first stage linear regression standard deviation
#'
#' @return returns a covariancne matrix


make.cov<- function(rho,tau0,tau1,y_sd,endog_sd){
  mat1 = matrix(c(1,rho,rho,y_sd^2),ncol = 2, byrow = T)
  tau_mat = as.matrix(cbind(unlist(tau0),unlist(tau1)))
  endog_mat = diag(length(endog_sd))*endog_sd^2
  Sig_err = rbind(cbind(mat1,t(tau_mat)),cbind(tau_mat,endog_mat))

  return(as.matrix(Sig_err))
}

#' Reconstitute a covariance matrix
#'
#' @description This takes a vector of upper diagonal terms and
#' reconstitutes it into the matrix you want in the end. The first
#'term is always 1 under the assumptions of the IV hurdle models.
#'If you did a cholesky decomposition it will un-do that tranformation
#'for you.
#'
#' @param vals a vector of covariance matrix elements
#' @param num the number of endogenous regressors
#' @param chol True/False. If true, will re-transform to original
#'
#' @return The vector renamed according to what it needs to be to work in the
#' likelihood functions


reconstitute.cov<-function(vals,num,chol=myChol){
  # de-cholesky-ify if you're supposed to
  if(chol == T){
    cov_vals = c(1,vals)
    empty = diag(2+num)
    empty[upper.tri(empty,diag = T)]<-cov_vals
    Sig_err = t(empty)%*%empty
    Sig_err = Sig_err/Sig_err[1,1]
  }
  else{
    print("Warning: no cholesky decomposition makes maximization less stable")
    cov_vals = c(1,vals)
    empty = diag(2+num)
    empty[upper.tri(empty,diag = T)]<-cov_vals
    tempty = t(empty)
    tempty[upper.tri(tempty,diag = T)]<-cov_vals
    Sig_err = tempty
  }
  return(Sig_err)
}

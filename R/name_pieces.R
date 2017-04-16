#' Name the elements that go into your regression function.
#'
#' @description There is likely
#' a cleaner way to do this, but it allows for more flexibility as to
#' how many of each parameter you're going to have.
#'
#' @param t a vector which should start with the covariates, then the
#' betas, then the gammas, then the pis
#' @param a1 the number of elements in the covariance matrix
#' @param a2 the number beta coefficients (2nd stage linear regression)
#' @param a3 the number endogenous regressors
#' @param a4 a list of the number of pi coefficients (1st stage regression)
#' for each endogenous regressor
#'
#' @return The vector renamed according to what it needs to be to work in the
#' likelihood functions

name.pieces<-function(t,a1 = len_cov,a2 = num_betas
                      ,a3 = num_endog,a4 = numpis){
  cov_start = t[1:a1]
  beta = t[(a1+1):(a1+a2)]
  gamma = t[(a1+a2+1):(a1+2*a2)]
  pis = list()
  k = a1+2*a2+1
  for(i in 1:a3){
    pis[[i]] = t[k:(k+a4[[i]]-1)]
    k = k+a4[[i]]
  }
  #pis = t[(a1+2*a2+1):length(t)]

  return(list(cov_start = cov_start
              ,beta = beta
              ,gamma = gamma
              ,pi = pis))
}

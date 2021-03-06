#' Lognormal IV likelihood
#'
#' @description Returns the likelihood for a lognormal IV. Requires certain
#' parameters to be attached/in the environment. Can't figure out how to
#' just pass them to the function without the function later trying to
#' optimize for them. This should essentially always be used in a wrapper
#' for that reason.
#'
#' @param t a vector of parameters to be maximized. Order should be:
#' covariance matrix elements (leaving out leading 1), then betas, then
#' gammas, then pis.
#'
#' @return returns the loglikelihood of the given parameter vector for the
#' dataset.


loglik_lgiv<-function(t){
  ############ preliminaries ##############
  # get everyone named
  if(length(t)!=len_cov + 2*num_betas + Reduce("+",numpis)){
    print(paste('num betas=',num_betas))
    print(paste('num pis=',numpis))
    print(paste('len cov=',len_cov))
    print(paste('t=',t))
    stop("Start vector should be len_cov + 2*num_betas + number of pis")
  }

  pieces = name.pieces(t)

  # reconstitute covariance matrix de-cholesky-ify if you're supposed to
  Sig_err = reconstitute.cov(vals=pieces$cov_start
                             ,num=num_endog
                             ,chol=myChol)

  # transform covariance matrix with betas and gammas
  if(is.null(names(pieces$gamma))){noname = T}
  else{noname = F}
  Sig = make.covTrans(Sig_err, num_endog, pieces$gamma, pieces$beta,noname = noname)
  if(all(is.na(Sig))){
    return(-Inf)
  }

  ################ get means ######################
  # get log of outcome and its censored values
  logy1 = log(outcome)
  censored = outcome<=0

  # get unconditional means of x2'sm
  n = length(pieces$pi)
  mu_x2 = matrix(c(rep(NA,n*length(outcome))),ncol = n)
  for(i in 1:n){
    mu_x2[,i]= ER_mat[[i]]%*%pieces$pi[[i]]
  }

  # get unconditional means of y's
  mu_y0 = y_mat%*%pieces$gamma
  mu_y1 = y_mat%*%pieces$beta

  # get conditional means
  k = dim(Sig)[1]
  j = num_endog
  sig2_x2 = as.matrix(Sig[3:k,3:k])

  #Parameters for y0star given x2
  if(cant.solve(sig2_x2)){return(-Inf)}
  mu_y0_x2 = mu_y0 + t(Sig[1,3:k]%*%solve(sig2_x2)%*%t(endog_mat-mu_x2))
  sig2_y0_x2 = Sig[1,1] - Sig[1,3:k]%*%solve(sig2_x2)%*%Sig[3:k,1]

  #Parameters for log(y1star) given x2
  mu_y1_x2 = mu_y1 + t(Sig[2,3:k]%*%solve(sig2_x2)%*%t(endog_mat-mu_x2))
  sig2_y1_x2 = Sig[2,2] - Sig[2,3:k]%*%solve(sig2_x2)%*%Sig[3:k,2]

  #Parameters for y0star given y1star and x2
  if(cant.solve(Sig[2:k,2:k])){return(-Inf)}
  mu_y0_y1x2 = mu_y0 + t(Sig[1,2:k,drop=FALSE]%*%solve(Sig[2:k,2:k])%*%t(cbind(logy1-mu_y1,endog_mat-mu_x2)))
  sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:k,drop=FALSE]%*%solve(Sig[2:k,2:k])%*%Sig[2:k,1,drop=FALSE]

  ############# calculate likelihood #############
  x2part = 0
  for(i in 1:j){
    temp = dnorm(endog_mat[,i],mean=mu_x2[,i],sd=sqrt(sig2_x2[i,i]),log = TRUE)
    x2part = x2part + temp
  }
  #When y1=0:
  ll0 = pnorm(0,mean=mu_y0_x2,sd=sqrt(sig2_y0_x2),log.p = TRUE) + x2part


  #When y1>0:
  if(sig2_y0_y1x2<0){return(-Inf)}

  ll1 = pnorm(0,mean=mu_y0_y1x2, sd=sqrt(sig2_y0_y1x2),log.p=TRUE, lower.tail = FALSE) +
    dnorm(logy1, mean=mu_y1_x2, sd=sqrt(sig2_y1_x2),log = TRUE) +
    x2part

  #Combine them, based on y1
  ll = ifelse(censored,ll0,ll1)

  #Return the negative log-likelihood

  -sum(ll, na.rm = T)
}



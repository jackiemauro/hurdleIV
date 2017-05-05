#' Cragg IV likelihood
#'
#' @description Returns the likelihood for a Cragg IV. Requires certain
#' parameters to be attached/in the environment. Can't figure out how to
#' just pass them to the function without the function later trying to
#' optimize for them. This should essentially always be used in a wrapper
#' for that reason. Note that this uses the first method of truncation for
#' the errors, described in cragg_errs.
#'
#' @param t a vector of parameters to be maximized. Order should be:
#' covariance matrix elements (leaving out leading 1), then betas, then
#' gammas, then pis.
#'
#' @return returns the loglikelihood of the given parameter vector for the
#' dataset.


loglik_craggiv<-function(t){
  ############ preliminaries ##############

  print("new round")
  require(mvtnorm, quietly = T)
  if(length(t)!=len_cov + 2*num_betas + Reduce("+",numpis)){
    stop("Start vector should be len_cov + 2*num_betas + number of pis")
  }

  # get everyone named
  pieces = name.pieces(t)

  # reconstitute covariance matrix de-cholesky-ify if you're supposed to
  Sig_err = reconstitute.cov(vals=pieces$cov_start
                             ,num=num_endog
                             ,chol=myChol)

  if(!all(Sig_err == t(Sig_err))){print("Sig_err not symmetric")}

  # transform covariance matrix with betas and gammas
  if(is.null(names(pieces$gamma))){noname = T}
  else{noname = F}

  Sig = make.covTrans(Sig_err, num_endog, pieces$gamma, pieces$beta,noname=noname)

  options(error = recover)
  tryCatch({
    sqrt(min(eigen(Sig)$value))
    },warning = function(w){
      print(paste("Sqrt warning: ",w))
      browser()
    }
  )

  ################ get means ######################
  # get censored values
  censored = outcome<=0

  # get unconditional means of x2's
  n = length(pieces$pi)
  mu_x2 = matrix(unlist(Map('%*%',ER_mat,pieces$pi)),ncol=n)

  # get unconditional means of y's
  mu_y0 = y_mat%*%pieces$gamma
  mu_y1 = y_mat%*%pieces$beta

  # get conditional means
  k = dim(Sig)[1]
  j = num_endog
  sig2_x2 = as.matrix(Sig[(k+1-j):k,(k+1-j):k])

  #Parameters for y0star given x2
  if(cant.solve(sig2_x2)){return(-Inf)}
  mu_y0_x2 = mu_y0 + t(Sig[1,(k+1-j):k]%*%solve(sig2_x2)%*%t(endog_mat-mu_x2))
  sig2_y0_x2 = Sig[1,1] - Sig[1,(k+1-j):k]%*%solve(sig2_x2)%*%Sig[(k+1-j):k,1]

  #Parameters for y1star given x2
  mu_y1_x2 = mu_y1 + t(Sig[2,(k+1-j):k]%*%solve(sig2_x2)%*%t(endog_mat-mu_x2))
  sig2_y1_x2 = Sig[2,2] - Sig[2,(k+1-j):k]%*%solve(sig2_x2)%*%Sig[(k+1-j):k,2]

  #Parameters for y0star given y1star and x2
  if(cant.solve(Sig[2:k,2:k])){return(-Inf)}
  mu_y0_y1x2 = mu_y0 + t(Sig[1,2:k,drop=FALSE]%*%solve(Sig[2:k,2:k])%*%t(cbind(outcome-mu_y1,endog_mat-mu_x2)))
  sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:k,drop=FALSE]%*%solve(Sig[2:k,2:k])%*%Sig[2:k,1,drop=FALSE]

  #Parameters for y0star and y1star given x2
  mu_y1y0_x2 = t(cbind(mu_y0,mu_y1)) + Sig[1:2,3:k]%*%solve(sig2_x2)%*%t(endog_mat-mu_x2)
  sig2_y1y0_x2 = Sig[1:2,1:2] - Sig[1:2,3:k]%*%solve(sig2_x2)%*%Sig[3:k,1:2]
  sig2_y1y0_x2[upper.tri(sig2_y1y0_x2)] <- sig2_y1y0_x2[lower.tri(sig2_y1y0_x2)]
  if(any(eigen(sig2_y1y0_x2)$value<0)){return(-Inf)}


  ############# calculate likelihood #############
  x2part = 0
  for(i in 1:j){
    temp = dnorm(endog_mat[,i],mean=mu_x2[,i],sd=sqrt(sig2_x2[i,i]),log = TRUE)
    x2part = x2part + temp
  }

  #probability y1>0 | x2
  #C = pnorm(0,mean = mu_y1_x2, sd = sqrt(sig2_y1_x2),log.p = T, lower.tail = FALSE)
  # should this just be P(y1>0)?
  C = pnorm(0, mean = mu_y1, sd = sqrt(Sig[2,2]), log.p = T, lower.tail = F)

  #ll0 integral
  l = c(0,-Inf); u = c(Inf,0)
  f <- function(x){pmvnorm(lower = l, upper = u, mean = x, sigma = sig2_y1y0_x2)}
  ll0_int <- tryCatch({
      apply(mu_y1y0_x2,2,f)
    }, error = function(e){
      print(paste("Error message from pmvnorm: ",e))
      print(sig2_y1y0_x2)
      return(Inf)
    }, warning = function(w){
      print(paste("Warning message from pmvnorm: ",w))
    }
    )

  #When y1=0:
  ll0 = x2part - C + ll0_int

  #When y1>0:
  ll1 = pnorm(0,mean=mu_y0_y1x2, sd=sqrt(sig2_y0_y1x2),log.p=TRUE, lower.tail = FALSE) +
    dnorm(outcome, mean=mu_y1_x2, sd=sqrt(sig2_y1_x2),log = TRUE) -
    C + x2part

  #Combine them, based on y1
  ll = ifelse(censored,ll0,ll1)

  #Return the negative log-likelihood
  print(paste("min eigen:",min(eigen(Sig)$value)))
  print(paste("likelihood: ",-sum(ll,na.rm =T)))
  -sum(ll, na.rm = T)
}




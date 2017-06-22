#' Cragg IV likelihood
#'
#' @description Returns the likelihood for a Cragg IV. Requires certain
#' parameters to be attached/in the environment. Note that this uses the
#' MG method for the errors, described in cragg_errs. This is following
#' the method worked out with Max.
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
  mu_x2 = c(matrix(unlist(Map('%*%',ER_mat,pieces$pi)),ncol=n))
  x2 = endog_mat$endog

  # get unconditional means of y's
  mu_y0 = c(y_mat%*%pieces$gamma)
  mu_y1 = c(y_mat%*%pieces$beta)

  # get conditional means
  k = dim(Sig)[1]
  if(k>3){print("Cragg IV can only handle one endogenous regressor at present");stop()}
  j = num_endog
  sig2_x2 = as.matrix(Sig[3:k,3:k])

  #Parameters for y0star given x2
  if(cant.solve(sig2_x2)){return(-Inf)}
  mu_y0_x2 = mu_y0 + t(Sig[1,3:k]%*%solve(sig2_x2)%*%t(x2-mu_x2))
  sig2_y0_x2 = Sig[1,1] - Sig[1,3:k]%*%solve(sig2_x2)%*%Sig[(k+1-j):k,1]

  #Parameters for y1star given x2
  mu_y1_x2 = mu_y1 + Sig[2,3]/Sig[3,3]*(x2-mu_x2)
  sig2_y1_x2 = Sig[2,2] - Sig[2,3]^2/Sig[3,3]

  #Parameters for y0star given y1star and x2
  mu_y0_y1x2 = c(mu_y0 + Sig[1,2:k,drop=FALSE]%*%solve(Sig[2:k,2:k])%*%rbind(y1-mu_y1,x2-mu_x2))
  sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:k,drop=FALSE]%*%solve(Sig[2:k,2:k])%*%Sig[2:k,1,drop=FALSE]

  #Parameters for y0star and y1star given x2
  mu_y1y0_x2 = t(cbind(mu_y0,mu_y1)) + Sig[1:2,3:k]%*%solve(sig2_x2)%*%t(x2-mu_x2)
  sig2_y1y0_x2 = Sig[1:2,1:2] - Sig[1:2,3:k]%*%solve(sig2_x2)%*%Sig[3:k,1:2]
  sig2_y1y0_x2[upper.tri(sig2_y1y0_x2)] <- sig2_y1y0_x2[lower.tri(sig2_y1y0_x2)]
  if(any(eigen(sig2_y1y0_x2)$value<0)){return(-Inf)}

  #Parameters for y1star given y0star and x2
  Sig_02 = Sig[-2,-2]
  sig2_y1_y0x2 = Sig[2,2] - Sig[2,c(1,3),drop=FALSE]%*%solve(Sig_02)%*%Sig[c(1,3),2,drop=FALSE]



  ############# calculate likelihood #############
  dat = data.frame(y1=y1,x2=x2
                   ,mu_y0_y1x2=mu_y0_y1x2
                   ,mu_y1_x2=mu_y1_x2,mu_y0_x2=mu_y0_x2
                   ,mu_x2=mu_x2,mu_y1=mu_y1,mu_y0=mu_y0)
  dat0 = dat[y1==0,]
  dat1 = dat[y1>0,]


  # ll0
  f0 <- function(dat){
    dat = as.list(dat)
    pnorm(0,mean=dat$mu_y0_x2,sd=sqrt(sig2_y0_x2),log.p = TRUE) +
      dnorm(dat$x2,mean=dat$mu_x2,sd=sqrt(sig2_x2),log = TRUE)
  }

  ll0 = sum( apply(dat0,1,f0))


  #ll1
  f1_max <- function(dat){
    prec=100
    dat = as.list(dat)
    counter_bad = 0

    mean_y1_y0x2 = function(y0,x2, mu_y1, mu_y0, mu_x2){mu_y1 + Sig[2,c(1,3),drop=FALSE]%*%solve(Sig_02)%*%rbind(y0-mu_y0,x2-mu_x2)}

    integrand <- function(u){
      top_func = function(u){
        mu_y1_y0x2 = mean_y1_y0x2(u,dat$x2, dat$mu_y1, dat$mu_y0, dat$mu_x2)
        dnorm.hack(dat$y1, mu_y1_y0x2, sqrt(sig2_y1_y0x2) ,log=TRUE) +
          dnorm.hack(u, dat$mu_y0_x2, sqrt(sig2_y0_x2),log=TRUE) + dnorm.hack(dat$x2, dat$mu_x2, sqrt(sig2_x2),log=TRUE)}

      bot_func = function(u){
        m = mean_y1_y0x2(u, dat$x2, dat$mu_y1, dat$mu_y0, dat$mu_x2)
        pnorm(mpfr(0,prec), m, rep(sqrt(sig2_y1_y0x2),length(m)), lower.tail=FALSE, log.p = TRUE)
      }

      top = top_func(u)
      bottom = bot_func(u)
      out = exp(top-bottom)
      out = as.numeric(out)
      out
    }

    int_result = integrate(integrand,0,Inf)

    log(int_result$value)
  }

  ll1 = sum( apply(dat1,1,f1_max) )

  # print(paste("neg likelihood=",-(ll0+ll1)))
  # print(paste('total time:',proc.time()[1] -tm0))

  -(ll0+ll1)
}




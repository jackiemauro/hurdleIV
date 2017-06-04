### Trying to do this MCMC style

loglik_dChoi = function(params){
  sig_v = params[1]
  sig_u = params[2]
  tau1 = params[3]
  tau0 = params[4]
  rho = params[5]
  beta11 = params[6]
  beta12 = params[7]
  beta2 = params[8]
  gamma11 = params[9]
  gamma12 = params[10]
  gamma2 = params[11]
  pi11 = params[12]
  pi12 = params[13]
  pi2 = params[14]

  censored = y1<=0

  pre = matrix( c(1, rho,    tau0,
                  0,     sig_u^2,    tau1,
                  0,  0, sig_v^2),
                ncol = 3, byrow = T)
  Sig_err = t(pre)%*%pre / ((t(pre)%*%pre)[1,1])
  if(min(eigen(Sig_err)$values)<=0){return(Inf)}
  if((sig_u<=0)|(sig_v<=0)){return(Inf)}

  #Transformation matrix from (eta, u, v) to (y0*, y1*, x2)
  A = rbind(
    c(1,0,gamma2),
    c(0,1,beta2),
    c(0,0,1)
  )

  #Covariance of (y0*, y1*, x2)
  Sig = A%*%Sig_err%*%t(A)
  if(min(eigen(Sig)$values)<=0){return(Inf)}

  #Means for (y0^*, y1*, x2)
  mu_y0 = gamma11 + gamma12*x1 + gamma2*(pi11 + pi12*x1 + pi2*z)
  mu_y1 = beta11 + beta12*x1 + beta2*(pi11 + pi12*x1 + pi2*z)
  mu_x2 = pi11 + x1*pi12 + z*pi2

  #Parameters for x2
  sig2_x2 = Sig[3,3]

  #Parameters for y0star given x2
  mu_y0_x2 = mu_y0 + Sig[1,3]/Sig[3,3]*(x2-mu_x2)
  sig2_y0_x2 = Sig[1,1] - Sig[1,3]^2/Sig[3,3]

  #Parameters for log(y1star) given x2
  mu_y1_x2 = mu_y1 + Sig[2,3]/Sig[3,3]*(x2-mu_x2)
  sig2_y1_x2 = Sig[2,2] - Sig[2,3]^2/Sig[3,3]

  #Parameters for y0star given y1star and x2
  mu_y0_y1x2 = mu_y0 + Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%rbind(y1-mu_y1,x2-mu_x2)
  sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%Sig[2:3,1,drop=FALSE]

  #Parameters for y0star and y1star given x2
  mu_y1y0_x2 = t(cbind(mu_y0,mu_y1)) + Sig[1:2,3]%*%solve(sig2_x2)%*%t(x2-mu_x2)
  sig2_y1y0_x2 = Sig[1:2,1:2] - Sig[1:2,3]%*%solve(sig2_x2)%*%Sig[3,1:2]
  sig2_y1y0_x2[upper.tri(sig2_y1y0_x2)] <- sig2_y1y0_x2[lower.tri(sig2_y1y0_x2)]
  if(any(eigen(sig2_y1y0_x2)$value<0)){return(-Inf)}

  # simulate denominator
  y0p = rnorm(1e6,mean=mu_y0,sd=sqrt(Sig[1,1]))
  y1p = rnorm(1e6,mean=mu_y1,sd=sqrt(Sig[2,2]))
  denom = mean(y1p>0&y0p>0)

  ###Calculate the contributions to the log likelihood
  #When y1=0:
  ll0 = pnorm(0,mean=mu_y0_x2,sd=sqrt(sig2_y0_x2),log.p = TRUE) +
    dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)

  #When y1>0: (not done)
  ll1 = pnorm(0,mean=mu_y0,sd=sqrt(Sig[1,1]),log.p = T, lower.tail = F)+
    pnorm(0,mean=mu_y0_y1x2,sd=sqrt(sig2_y0_y1x2),log.p = T,lower.tail = F)+
    dnorm(y1, mean=mu_y1_x2, sd=sqrt(sig2_y1_x2),log = TRUE) +
    dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)-
    log(denom)

  #Combine them, based on y1
  ll = ifelse(censored,ll0,ll1)
  print(-sum(ll))

  #Return the negative log-likelihood
  -sum(ll)
}

### This function returns the negative log-likelihood (so that it can be minimized)
### The math follows the derivation in the pdf I sent over
### There is not much error checking on the parameters, so be careful.
### The function does check for a PSD (eta,u,v) covariance matrix, and returns Inf if it fails.  This works well with optim.
### Now just adding a covariance parameter instead of assuming independence u,eta


### this function was created by Max, it is simpler so easier for testing.
loglik_lgnorm_rho = function(params){
  sig_v = params[1]
  sig_u = params[2]
  tau1 = params[3]
  tau0 = params[4]
  beta1 = params[5]
  beta2 = params[6]
  gamma1 = params[7]
  gamma2 = params[8]
  pi1 = params[9]
  pi2 = params[10]
  rho = params[11]

  censored = y1<=0
  #We will just work with the log of y1, so that we don't need to care about lognormal distributions.
  logy1 = log(y1)

  #Check for invalid parameters
  Sig_err= matrix( c(1, rho,    tau0,
                     rho,     sig_u^2,    tau1,
                     tau0,  tau1, sig_v^2),
                   ncol = 3, byrow = T)

  if(min(eigen(Sig_err)$values)<=0){return(Inf)}
  if((sig_u<=0)|(sig_v<=0)){return(Inf)}

  #Transformation matrix from (eta, u, v) to (y0*, log y1*, x2) (without the mean)
  A = rbind(
    c(1,0,gamma2),
    c(0,1,beta2),
    c(0,0,1)
  )
  #Covariance of (y0*, log y1*, x2)
  Sig = A%*%Sig_err%*%t(A)
  x2 = as.numeric(x2)#I was having some issues with the format of x2 before, maybe this doesn't matter anymore.

  #### WARNING: I WILL BE SLOPPY FROM HERE ON BY DROPPING THE LOG IN MY log y1* NOTATION (and the star) ####
  #### It's just too much to carry around

  #Means for (y0^*, log y1*, x2).
  mu_y0 = (gamma1+pi1*gamma2)*x1 + gamma2*pi2*z
  mu_y1 = (beta1+pi1*beta2)*x1+pi2*beta2*z
  mu_x2 = x1*pi1 + z*pi2

  #### NOTATION: For mu and sig2 variables, I will separate the variables from the conditioning variables by underscores.
  #### That is, (mu/sig2)_(Variable)_(Conditioning variables)

  #Parameters for x2
  sig2_x2 = Sig[3,3]

  #Parameters for y0star given x2
  mu_y0_x2 = mu_y0 + Sig[1,3]/Sig[3,3]*(x2-mu_x2)
  sig2_y0_x2 = Sig[1,1] - Sig[1,3]^2/Sig[3,3]

  #Parameters for log(y1star) given x2
  mu_y1_x2 = mu_y1 + Sig[2,3]/Sig[3,3]*(x2-mu_x2)
  sig2_y1_x2 = Sig[2,2] - Sig[2,3]^2/Sig[3,3]

  #Parameters for y0star given y1star and x2
  mu_y0_y1x2 = mu_y0 + Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%rbind(logy1-mu_y1,x2-mu_x2)
  sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%Sig[2:3,1,drop=FALSE]

  ###Calculate the contributions to the log likelihood.

  #When y1=0:
  ll0 = pnorm(0,mean=mu_y0_x2,sd=sqrt(sig2_y0_x2),log.p = TRUE) + dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)

  #When y1>0:
  ll1 = pnorm(0,mean=mu_y0_y1x2, sd=sqrt(sig2_y0_y1x2),log.p=TRUE, lower.tail = FALSE) +
    dnorm(logy1, mean=mu_y1_x2, sd=sqrt(sig2_y1_x2),log = TRUE) +
    dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)

  #Combine them, based on y1
  ll = ifelse(censored,ll0,ll1)

  #Return the negative log-likelihood, since I will be using optim (which minimizes instead of maximizing)
  -sum(ll)
}

# ### this function uses dmvnorm. Putting it this way for MCMC because similar tests performed better
# loglik_lgnorm_dmv = function(params){
#   sig_v = params[1]
#   sig_u = params[2]
#   tau1 = params[3]
#   tau0 = params[4]
#   beta1 = params[5]
#   beta2 = params[6]
#   gamma1 = params[7]
#   gamma2 = params[8]
#   pi1 = params[9]
#   pi2 = params[10]
#   rho = params[11]
#
#   censored = y1<=0
#   #We will just work with the log of y1, so that we don't need to care about lognormal distributions.
#   logy1 = log(y1)
#
#   #Check for invalid parameters
#   Sig_err= matrix( c(1, rho,    tau0,
#                      rho,     sig_u^2,    tau1,
#                      tau0,  tau1, sig_v^2),
#                    ncol = 3, byrow = T)
#   if(min(eigen(Sig_err)$values)<=0){return(Inf)}
#   if((sig_u<=0)|(sig_v<=0)){return(Inf)}
#
#   #Transformation matrix from (eta, u, v) to (y0*, log y1*, x2) (without the mean)
#   A = rbind(
#     c(1,0,gamma2),
#     c(0,1,beta2),
#     c(0,0,1)
#   )
#   #Covariance of (y0*, log y1*, x2)
#   Sig = A%*%Sig_err%*%t(A)
#   x2 = as.numeric(x2)#I was having some issues with the format of x2 before, maybe this doesn't matter anymore.
#
#   #### WARNING: I WILL BE SLOPPY FROM HERE ON BY DROPPING THE LOG IN MY log y1* NOTATION (and the star) ####
#   #### It's just too much to carry around
#
#   #Means for (y0^*, log y1*, x2).
#   mu_y0 = (gamma1+pi1*gamma2)*x1 + gamma2*pi2*z
#   mu_y1 = (beta1+pi1*beta2)*x1+pi2*beta2*z
#   mu_x2 = x1*pi1 + z*pi2
#
#   #### NOTATION: For mu and sig2 variables, I will separate the variables from the conditioning variables by underscores.
#   #### That is, (mu/sig2)_(Variable)_(Conditioning variables)
#
#   #Parameters for x2
#   sig2_x2 = Sig[3,3]
#
#   #Parameters for y0star given x2
#   mu_y0_x2 = mu_y0 + Sig[1,3]/Sig[3,3]*(x2-mu_x2)
#   sig2_y0_x2 = Sig[1,1] - Sig[1,3]^2/Sig[3,3]
#
#   #Parameters for log(y1star) given x2
#   mu_y1_x2 = mu_y1 + Sig[2,3]/Sig[3,3]*(x2-mu_x2)
#   sig2_y1_x2 = Sig[2,2] - Sig[2,3]^2/Sig[3,3]
#
#   #Parameters for y0star given y1star and x2
#   mu_y0_y1x2 = mu_y0 + Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%rbind(logy1-mu_y1,x2-mu_x2)
#   sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%Sig[2:3,1,drop=FALSE]
#
#   ###Calculate the contributions to the log likelihood.
#
#   #When y1=0:
#   ll0 = pnorm(0,mean=mu_y0_x2,sd=sqrt(sig2_y0_x2),log.p = TRUE) + dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)
#
#   #When y1>0:
#   mat = data.frame(x=c(x2),y=c(logy1),predx = c(mu_x2), predy = c(mu_y1))
#   singlelikelihoods = apply(mat,1,function(x) dmvnorm(x[1:2],mean = x[3:4],sigma = Sig_err[-1,-1],log = T))
#   ll1 = pnorm(0,mean=mu_y0_y1x2, sd=sqrt(sig2_y0_y1x2),log.p=TRUE, lower.tail = FALSE) +
#     sum(singlelikelihoods)
#
#   #Combine them, based on y1
#   ll = ifelse(censored,ll0,ll1)
#
#   #Return the negative log-likelihood, since I will be using optim (which minimizes instead of maximizing)
#   -sum(ll)
# }
#

# another attempt at mcmc

#### laplace demon http://www.sumsar.net/blog/2013/06/three-ways-to-run-bayesian-models-in-r/####
rm(list = ls())
library(LaplacesDemon)

y = rnorm(1000,4,2)
# The model specification
model <- function(parm, data) {
  mu <- parm[1]
  sigma <- exp(parm[2])
  log_lik <- sum(dnorm(data$y, mu, sigma, log = T))
  post_lik <- log_lik + dnorm(mu, 0, 100, log = T) + dlnorm(sigma, 0, 4, log = T)
  # This list is returned and has to follow a the format specified by
  # LaplacesDemon.
  list(LP = post_lik, Dev = -2 * log_lik, Monitor = c(post_lik, sigma), yhat = NA,
       parm = parm)
}

# Running the model
data_list <- list(N = length(y), y = y, mon.names = c("post_lik", "sigma"),
                  parm.names = c("mu", "log.sigma"))
mcmc_samples <- LaplacesDemon(Model = model, Data = data_list, Iterations = 30000,
                              Algorithm = "HARM", Thinning = 1)

plot(mcmc_samples, BurnIn = 10000, data_list)

# basic IV
rm(list = ls())
require(mvtnorm); require(MASS)
N = 10000
z <- rnorm(N,2,1)
sig2y1<-2; sig2x<-3; t0<-.2
Pi1<-1; Pi2<-3; B1<-2; B2<-1
errs = mvrnorm(n = N, mu = c(0,0), Sigma = matrix(c(sig2y1,t0,t0,sig2x),ncol = 2, byrow = T))
x2 <- Pi1 + Pi2*z + errs[,2]
y1 <- B1 + B2*x2 + errs[,1]

model <- function(parm, data) {
  b1 <- parm[1]; b2 <- parm[2]; p1 <- parm[3]; p2 <- parm[4]
  sigy1 <- exp(parm[5]); sigx <- exp(parm[6]); t0 <- parm[7]

  mu_x <- p1 + p2*data$z
  mu_y <- b1 + b2*mu_x
  Sig = matrix(c(sigy1,t0,t0,sigx),ncol = 2, byrow = T)

  #Transformation matrix from (eta, u, v) to (y0*, log y1*, x2) (without the mean)
  A = rbind(c(1,b2),c(0,1))
  #Covariance of (y0*, log y1*, x2)
  Sig2 = A%*%Sig%*%t(A)

  mu_y_x = mu_y + (Sig2[1,2]/Sig2[2,2])*(data$x-mu_x)
  sig2_y_x = Sig2[1,1] - (Sig2[1,2]^2)/Sig2[2,2]

  if(sig2_y_x>0){
    log_lik <- sum(dnorm(data$y,mu_y_x,sqrt(sig2_y_x),log = T) + dnorm(data$x,mu_x,sqrt(Sig2[2,2]),log =T))
  }
  else{
    mat = data.frame(y=data$y,x=data$x,predy=mu_y,predx = mu_x)
    singlelikelihoods = apply(mat,1,function(x) dmvnorm(x[1:2],mean = x[3:4],sigma = Sig,log = T))
    log_lik <- sum(singlelikelihoods)
  }

  log_prior <- sum(dnorm(c(b1,b2,p1,p2,t0),0,100,log = T))+dlnorm(sigy1, 0, 4, log = T)+dlnorm(sigx, 0, 4, log = T)

  post_lik <- log_lik + log_prior

  # This list is returned and has to follow a the format specified by
  # LaplacesDemon.
  list(LP = post_lik, Dev = -2 * log_lik, Monitor = c(post_lik, sigx,sigy1), yhat = NA,parm = parm)
}

# Running the model
data_list <- list(N = length(y1), y = y1,x = x2, z = z , mon.names = c("post_lik", "sigx","sigy1"),
                  parm.names = c("B1","B2","Pi1","Pi2","log.sigy","log.sigx","tau0"))
mcmc_samples <- LaplacesDemon(Model = model, Data = data_list, Iterations = 50000,
                              Algorithm = "HARM", Thinning = 1)

plot(mcmc_samples, BurnIn = 10000, data_list)


# lognormal IV
rm(list = ls())
N = 10000
z <- rnorm(N,2,1); x1 <- rnorm(N,-1,3)
sig2y1<-2; sig2y0<-1; rho <- .1; t0<-.2; t1 <- .3; sig2x<-3
Pi1<-1; Pi2<-3; B1<-2; B2<-1;G1<- -1; G2 <- 2
errs = mvrnorm(n = N, mu = c(0,0,0), Sigma = matrix(c(sig2y0,rho,t0,rho,sig2y1,t1,t0,t1,sig2x),ncol = 3, byrow = T))
x2 <- Pi1*x1 + Pi2*z + errs[,3]
y0Star <- G1*x1 + G2*z + errs[,2]
lny1Star <- B1*x1 + B2*x2 + errs[,1]
y1 <- as.numeric(y0Star>0)*exp(lny1Star)

model <- function(parm, data) {
  sigx <- exp(parm[1]); sigy1 <- exp(parm[2]); t1 <- parm[3]; t0 <- parm[4]
  b1 <- parm[5]; b2 <- parm[6]; g1 <- parm[7]; g2 <- parm[8]; p1<-parm[9];p2<-parm[10];rho <- parm[11]

  log_lik <- loglik_lgnorm_rho(parm)

  log_prior <- sum(dnorm(c(b1,b2,g1,g2,p1,p2,t1,t0,rho),0,100,log = T))+dlnorm(sigy1, 0, 4, log = T)+dlnorm(sigx, 0, 4, log = T)

  post_lik <- log_lik + log_prior
  print(post_lik)

  # This list is returned and has to follow a the format specified by
  # LaplacesDemon.
  list(LP = post_lik, Dev = -2 * log_lik, Monitor = c(post_lik, sigx,sigy1), yhat = NA,parm = parm)
}

# Running the model -- infinite eigen(sig_err values)
data_list <- list(N = length(y1), y = y1,x = x2, x1=x1,z = z , mon.names = c("post_lik", "sigx","sigy1"),
                  parm.names = c("log.sigx","log.sigy","tau1","tau0","B1","B2","G1","G2","Pi1","Pi2","rho"))
mcmc_samples <- LaplacesDemon(Model = model, Data = data_list, Iterations = 50000,
                              Algorithm = "HARM", Thinning = 1,
                              Initial.Values = c(1,1,0,0,B1,B2,G1,G2,Pi1,Pi2,rho))

plot(mcmc_samples, BurnIn = 10000, data_list)


# Running the model
data_list <- list(N = length(y1), y1 = y1, x2 = x2, z = z, mon.names = c("post_lik", "sig_v","sig_u"),
                  parm.names = c("log.sigv","log.sigu", "t1", "t0","b1","b2","g1","g2", "p1", "p2"))
mcmc_samples <- LaplacesDemon(Model = model, Data = data_list, Iterations = 30000,
                              Algorithm = "HARM", Thinning = 1)

plot(mcmc_samples, BurnIn = 10000, data_list)

# ##### this is from: https://nicercode.github.io/guides/mcmc/ ####
# just to get mean
m = 0
f <- function(x) dnorm(x,m,1)
q <- function(x) rnorm(1,x,1)
step <- function(x,f,q){
  xp = q(x) #propose a new x
  alpha = min(1,f(xp)/f(x)) #calculate alpha
  if(runif(1)<alpha){x<-xp}
  x
}

run <- function(x, f, q, nsteps) {
  res <- matrix(NA, nsteps, length(x))
  for (i in seq_len(nsteps))
    res[i,] <- x <- step(x, f, q)
  drop(res)
}

res <- run(-10, f, q, 1000)

layout(matrix(c(1, 2), 1, 2), widths=c(4, 1))
par(mar=c(4.1, .5, .5, .5), oma=c(0, 4.1, 0, 0))
plot(res, type="s", xpd=NA, ylab="Parameter", xlab="Sample", las=1)
usr <- par("usr")
xx <- seq(usr[3], usr[4], length=301)
plot(f(xx), xx, type="l", yaxs="i", axes=FALSE, xlab="")

#### looking at these guys' example: https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/####
rm(list = ls())
trueA <- 5
trueB <- 0
trueSd <- 10
sampleSize <- 31

# create independent x-values
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
# create dependent values according to ax + b + N(0,sd)
y <-  trueA * x + trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)

likelihood <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]

  pred = a*x + b
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)
}

prior <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  aprior = dunif(a, min=0, max=10, log = T)
  bprior = dnorm(b, sd = 5, log = T)
  sdprior = dunif(sd, min=0, max=30, log = T)
  return(aprior+bprior+sdprior)
}

posterior <- function(param){
  return (likelihood(param) + prior(param)) #being bayesian
  #return(likelihood(param))#what if i don't be bayesian?
}

proposalfunction <- function(param){
  return(rnorm(3,mean = param, sd= c(0.1,0.5,0.3)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])

    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

startvalue = c(4,0,10)
chain1 = run_metropolis_MCMC(startvalue, 10000) #bayesian
chain2 = run_metropolis_MCMC(startvalue, 10000) #not bayesian (change posterior)

burnIn = 5000
acceptance = 1-mean(duplicated(chain1[-(1:burnIn),]))

names = c('a','b','sig'); true = c(trueA,trueB,trueSd)
for(i in 1:dim(chain1)[2]){
  hist(chain1[-(1:burnIn),i],nclass=30, main=paste("Posterior of ",names[i]), xlab="True value = red line" )
  abline(v = mean(chain1[-(1:burnIn),i]));abline(v = true[i], col="red" )
}

#### basic IV #####
# even though dmvnorm and dnorm using conditionals give the same probabilities
# asking the algorithm to estimate the mean of x/y within the likelihood makes it
# perform less well

rm(list = ls())
require(mvtnorm)
N = 10000
z <- rnorm(N,2,1)
sig2y1<-2; sig2x<-3; t0<-.2
Pi1<-1; Pi2<-3; B1<-2; B2<-1
errs = mvrnorm(n = N, mu = c(0,0), Sigma = matrix(c(sig2y1,t0,t0,sig2x),ncol = 2, byrow = T))
x <- Pi1 + Pi2*z + errs[,2]
y <- B1 + B2*x + errs[,1]

likelihood <- function(param){
  s1 = param[1];s2 = param[2]; t = param[3]; p1 = param[4]; p2 = param[5]; b1 = param[6]; b2 = param[7]

  xpred = p1 + p2*z
  ypred = b1 + b2*xpred
  preds = matrix(c(ypred,xpred),ncol = 2, byrow = F)
  obs = matrix(c(y,x),ncol = 2, byrow = F)
  Sig = matrix(c(s1,t,t,s2),ncol = 2, byrow = T)

  # These two will give the same likelihoods, except, I think, when sig2_y1_x2<0
  # Sig_err = matrix(c(s1,t,t,s2),ncol = 2, byrow = T)
  # A = matrix(c(1,b1,0,1),ncol = 2, byrow = T)
  # Sig = A%*%Sig_err%*%t(A)
  #
  # mu_y1_x2 = ypred + Sig[1,2]/Sig[2,2]*(x-xpred)
  # sig2_y1_x2 = Sig[1,1] - Sig[1,2]^2/Sig[2,2]
  #
  # ylik = sum(dnorm(y,mean = mu_y1_x2,sd = sqrt(sig2_y1_x2),log = T))
  # xlik = sum(dnorm(x,mean = xpred,sd = sqrt(Sig[2,2]),log = T))
  # sumll1 = xlik + ylik
  #
  mat = data.frame(y=y,x=x,predy=ypred,predx = xpred)
  singlelikelihoods = apply(mat,1,function(x) dmvnorm(x[1:2],mean = x[3:4],sigma = Sig,log = T))
  sumll = sum(singlelikelihoods)
  #
  # print(sumll)
  # print(sumll1)

  return(sumll)
}

likelihood <- function(param){
  # this one has to guess at the true value of x in order to get the conditional
  # mean of y--faster but appears not to be as accurate.
  s1 = param[1];s2 = param[2]; t = param[3]; p1 = param[4]; p2 = param[5]; b1 = param[6]; b2 = param[7]

  xpred = p1 + p2*z
  ypred = b1 + b2*xpred
  preds = matrix(c(xpred,ypred),ncol = 2, byrow = F)
  obs = matrix(c(x,y),ncol = 2, byrow = F)
  #Sig = matrix(c(s1,t,t,s2),ncol = 2, byrow = T)
  Sig_err = matrix(c(s1,t,t,s2),ncol = 2, byrow = T)
  A = matrix(c(1,b1,0,1),ncol = 2, byrow = T)
  Sig = A%*%Sig_err%*%t(A)

  mu_y1_x2 = ypred + Sig[1,2]/Sig[2,2]*(x-xpred)
  sig2_y1_x2 = Sig[1,1] - Sig[1,2]^2/Sig[2,2]
  if(sig2_y1_x2<=0 | Sig[2,2]<=0){return(-Inf)}

  ylik = sum(dnorm(y,mean = mu_y1_x2,sd = sqrt(sig2_y1_x2),log = T))
  xlik = sum(dnorm(x,mean = xpred,sd = sqrt(Sig[2,2]),log = T))
  sumll = xlik + ylik
  return(sumll)
}


prior <- function(param){
  s1 = param[1];s2 = param[2]; t = param[3]; p1 = param[4]; p2 = param[5]; b1 = param[6]; b2 = param[7]

  s1prior = dunif(s1, min=0, max=10, log = T)
  s2prior = dunif(s2, min=0, max=10, log = T)
  tprior = dnorm(t, sd = 5, log = T)
  p1prior = dnorm(p1, sd = 5, log = T)
  p2prior = dnorm(p1, sd = 5, log = T)
  b1prior = dnorm(p1, sd = 5, log = T)
  b2prior = dnorm(p1, sd = 5, log = T)
  return(s1prior + s2prior + tprior + p1prior + p2prior + b1prior + b2prior)
}

posterior <- function(param){
  return (likelihood(param) + prior(param)) #being bayesian
  #return(likelihood(param))#what if i don't be bayesian?
}

proposalfunction <- function(param){
  return(rnorm(7,mean = param, sd= c(0.1,0.5,0.3,.5,.5,.5,.5)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,length(startvalue)))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])

    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (!is.na(probab) & (runif(1) < probab)){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

names = c('sig2y1','sig2x','tau','pi1','pi2','b1','b2')
true = c(sig2y1,sig2x,t0,Pi1,Pi2,B1,B2)

startvalue = c(rep(1,length(true)))
chain1 = run_metropolis_MCMC(startvalue, 1000)
chain2 = run_metropolis_MCMC(startvalue, 1000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain1[-(1:burnIn),]))


for(i in 1:dim(chain1)[2]){
  hist(chain1[,i],nclass=30, main=paste("Posterior of ",names[i]), xlab="True value = red line" )
  abline(v = mean(chain1[,i]));abline(v = true[i], col="red" )
}

xmins1 = apply(chain1,2,min);xmins2 = apply(chain2,2,min);xmins = apply(cbind(xmins1,xmins2),1,min)
xmaxs1 = apply(chain1,2,max);xmaxs2 = apply(chain2,2,max);xmaxs = apply(cbind(xmaxs1,xmaxs2),1,max)

par(mfrow = c(1,2))
for(i in 1:dim(chain2)[2]){
  hist(chain1[,i],nclass=30, main=paste("Posterior of ",names[i]," chain 1"), xlim = c(xmins[i],xmaxs[i]),xlab="True value = red line" )
  abline(v = mean(chain1[,i]));abline(v = true[i], col="red" )

  hist(chain2[,i],nclass=30, main=paste("Posterior of ",names[i], "chain 2"), xlim = c(xmins[i],xmaxs[i]), xlab="True value = red line" )
  abline(v = mean(chain2[,i]));abline(v = true[i], col="red" )
}
par(mfrow = c(1,1))

#### with fake loglik lgnorm ####
rm(list = ls())
require(mvtnorm)
N = 10000
z <- rnorm(N,2,1)
sig2y1<-2; sig2y0<-1; rho <- .1; t0<-.2; t1 <- .3; sig2x<-3
Pi1<-1; Pi2<-3; B1<-2; B2<-1;G1<- -1; G2 <- 2
errs = mvrnorm(n = N, mu = c(0,0,0), Sigma = matrix(c(sig2y0,rho,t0,rho,sig2y1,t1,t0,t1,sig2x),ncol = 3, byrow = T))
x <- Pi1 + Pi2*z + errs[,3]
y0Star <- G1 + G2*z + errs[,2]
y1Star <- B1 + B2*x + errs[,1]
y1 <- as.numeric(y0Star>0)*y1Star

likelihood <- function(param){ #not done
  sy0 = param[1]; sy1 = param[2]; sx2 = param[3];
  r = param[4]; t0 = param[5]; t1 = param[6]
  p1 = param[7]; p2 = param[8];g1 = param[9];g2 = param[10]; b1 = param[11]; b2 = param[12]

  xpred = p1 + p2*z
  y0pred = g1 + g2*xpred
  y1pred = b1 + b2*xpred
  preds = matrix(c(ypred,xpred),ncol = 2, byrow = F)
  obs = matrix(c(y,x),ncol = 2, byrow = F)
  Sig = matrix(c(sy0,r,t0,r,sy1,t1,t0,t1,sx2),ncol = 3, byrow = T)

  mat = data.frame(y=y,x=x,predy=ypred,predx = xpred)
  singlelikelihoods = apply(mat,1,function(x) dmvnorm(x[1:2],mean = x[3:4],sigma = Sig,log = T))
  sumll = sum(singlelikelihoods)

  return(sumll)
}


prior <- function(param){
  s1 = param[1];s2 = param[2]; t = param[3]; p1 = param[4]; p2 = param[5]; b1 = param[6]; b2 = param[7]

  s1prior = dunif(s1, min=0, max=10, log = T)
  s2prior = dunif(s2, min=0, max=10, log = T)
  tprior = dnorm(t, sd = 5, log = T)
  p1prior = dnorm(p1, sd = 5, log = T)
  p2prior = dnorm(p1, sd = 5, log = T)
  b1prior = dnorm(p1, sd = 5, log = T)
  b2prior = dnorm(p1, sd = 5, log = T)
  return(s1prior + s2prior + tprior + p1prior + p2prior + b1prior + b2prior)
}

posterior <- function(param){
  return (likelihood(param) + prior(param)) #being bayesian
  #return(likelihood(param))#what if i don't be bayesian?
}

proposalfunction <- function(param){
  return(rnorm(7,mean = param, sd= c(0.1,0.5,0.3,.5,.5,.5,.5)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,length(startvalue)))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])

    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (!is.na(probab) & (runif(1) < probab)){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

names = c('sig2y1','sig2x','tau','pi1','pi2','b1','b2')
true = c(sig2y1,sig2x,t0,Pi1,Pi2,B1,B2)

startvalue = c(rep(1,length(true)))
chain1 = run_metropolis_MCMC(startvalue, 1000)
chain2 = run_metropolis_MCMC(startvalue, 1000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain1[-(1:burnIn),]))


for(i in 1:dim(chain1)[2]){
  hist(chain1[,i],nclass=30, main=paste("Posterior of ",names[i]), xlab="True value = red line" )
  abline(v = mean(chain1[,i]));abline(v = true[i], col="red" )
}

#### with loglik_lgnorm ####
rm(list = ls())
detach(dat)
ps = c(0,-1,3); gs = c(0,.8,.07); bs = c(0,.06,.02); sigv = 5; sigu = 2;r=.2;t0=.3;t1=.1
true = c(sigv,sigu,t0,t1,bs,gs,ps,r)
dat = hurdle.IV.sim(pi=ps,gamma = gs,beta = bs,endog_sd = sigv,y_sd = sigu,rho =r,tau0 = t0,tau1=t1)
dat$y1 <-dat$y
dat$x2<-dat$endog;dat$x1<-dat$exog1;dat$z<-dat$inst1
attach(dat)

likelihood<-function(param) {-loglik_lgnorm_rho(param)}

prior <- function(params){
  #need to get this to be wishart
  invsds = dnorm(params[1:2],sd = 100,log = T)
  sds = 1/invsds
  lasts = Reduce("+",lapply(params[3:length(params)],function(x) dnorm(x,100,log=T)))
  lasts + sum(sds)
}

posterior <- function(param){
  return (likelihood(param) + prior(param)) #being bayesian
  #return(likelihood(param))#what if i don't be bayesian?
}

proposalfunction <- function(param){
  return(rnorm(length(param),mean = param, sd= rep(1,length(param))))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,length(startvalue)))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])

    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if ((!is.na(probab)) & runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

chain1 = run_metropolis_MCMC(true, 1e8) #bayesian from true
#chain1 = run_metropolis_MCMC(c(rep(1,14)),10000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain1[-(1:burnIn),]))

true = c(sigv,sigu,t0,t1,bs,gs,ps,r)
nm_pars = c('sigv','sigu','tau0','tau1','beta1','beta2','beta3','gamma1','gamma2','gamma3','pi1','pi2','pi3','rho')

for(i in 1:length(true)){
  hist(chain1[-(1:burnIn),i],nclass=30, main=paste("Posterior of",nm_pars[i])
       , xlab="True value = red line" )
  abline(v = mean(chain1[-(1:burnIn),i]));abline(v = true[i], col="red" )
}

# ### This function returns the negative log-likelihood (so that it can be minimized)
# ### The math follows Max G'Sell's derivation.
# 
# y1 = 0.04016536
# x2 = 5.863038
# x1 = 4.138829
# z = 0.650185
# 
# beta11 = 0.636857
# beta12 = -0.2093
# beta2 = -0.008241487
# 
# gamma11 = 0.343783236
# gamma12 = 0.267480127
# gamma2 = 0.119696815
# 
# 
# pi11 = -  1.067030153
# pi12 = 1.567173444
# pi2 = 1.532754633
# 
# sig_v = 1.714909010
# sig_u = 0.072066543
# tau1 = 0.054132764
# tau0 = 0.415410701
# rho = 0.804469355
# 
# params = c(sig_v,sig_u,tau1,tau0,rho,sig_v,sig_u,tau1,tau0,rho,beta11,beta12,beta2,gamma11,gamma12,gamma2,pi11,pi12,pi2)
#
#Not working yet
#library(parallel)
#clust = makeCluster(getOption("cl.cores", 2))


dnorm.hack = function (x, mean = 0, sd = 1, log = FALSE) 
{
  # if (is.numeric(x) && is.numeric(mean) && is.numeric(sd)) 
  #   stats__dnorm(x, mean, sd, log = log)
  # else if ((x.mp <- is(x, "mpfr")) || is(mean, "mpfr") || (s.mp <- is(sd, "mpfr"))) {
  #   prec <- pmax(53, getPrec(x), getPrec(mean), getPrec(sd))
  #   if (!x.mp) 
  #     x <- mpfr(x, prec)

  #   if (!s.mp) 
  #     sd <- mpfr(sd, prec)
  #   if (log) 
  #     -(log(sd) + (log(twopi) + x * x)/2)
  #   else exp(-x^2/2)/(sd * sqrt(twopi))
  # }
  # else stop("invalid arguments (x,mean,sd)")
  prec <- pmax(53, getPrec(x), getPrec(mean), getPrec(sd))
  n = max(length(x), length(mean),length(sd))
  x <- mpfr(x, prec)
  mean = mpfr(mean,prec)
  sd = mpfr(sd, prec)

  x <- (x - mean)/rep(sd,n)
  twopi <- 2 * Const("pi", prec)
  #browser()
  if (log) 
    -(rep(log(sd),n) + (log(twopi) + x * x)/2)
  else exp(-x^2/2)/(sd * sqrt(twopi))
}
  



loglik_MG = function(params){
  require(Rmpfr)
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
  tm0 = proc.time()[1]

  
  
  censored = y1<=0

  pre = matrix( c(1, rho,    tau0,
                  0,     sig_u^2,    tau1,
                  0,  0, sig_v^2),
                ncol = 3, byrow = T)
  Sig_err = t(pre)%*%pre / ((t(pre)%*%pre)[1,1])
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
  if(min(eigen(Sig)$values)<=0){return(Inf)}

  
  
  
  #Means for (y0^*, y1*, x2).
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
  mu_y0_y1x2 = c(mu_y0 + Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%rbind(y1-mu_y1,x2-mu_x2))
  sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%Sig[2:3,1,drop=FALSE]

  #Parameters for y0star and y1star given x2
  #mu_y1y0_x2 = c(t(t(cbind(mu_y0,mu_y1)) + Sig[1:2,3]%*%solve(sig2_x2)%*%t(x2-mu_x2)))
  sig2_y1y0_x2 = Sig[1:2,1:2] - Sig[1:2,3]%*%solve(sig2_x2)%*%Sig[3,1:2]
  sig2_y1y0_x2[upper.tri(sig2_y1y0_x2)] <- sig2_y1y0_x2[lower.tri(sig2_y1y0_x2)]
  if(any(eigen(sig2_y1y0_x2)$value<0)){return(-Inf)}

  #Parameters for y1star and y0star given x2
  Sig_02 = Sig[-2,-2]
  #mu_y1_y0x2 = mu_y1 + Sig[2,c(1,3),drop=FALSE]%*%solve(Sig_02)%*%rbind(y0-mu_y0,x2-mu_x2)
  sig2_y1_y0x2 = Sig[2,2] - Sig[2,c(1,3),drop=FALSE]%*%solve(Sig_02)%*%Sig[c(1,3),2,drop=FALSE]


  ###Calculate the contributions to the log likelihood.
  dat = data.frame(y1=y1,x2=x2
                   #,mu_y1y0_x2=mu_y1y0_x2
                   ,mu_y0_y1x2=mu_y0_y1x2
                   ,mu_y1_x2=mu_y1_x2,mu_y0_x2=mu_y0_x2
                   ,mu_x2=mu_x2,mu_y1=mu_y1,mu_y0=mu_y0)#, mu_y1_y0x2=mu_y1_y0x2)
  dat0 = dat[y1==0,]
  dat1 = dat[y1>0,]


  # ll0
  f0 <- function(dat){
    dat = as.list(dat)
    pnorm(0,mean=dat$mu_y0_x2,sd=sqrt(sig2_y0_x2),log.p = TRUE) +
      dnorm(dat$x2,mean=dat$mu_x2,sd=sqrt(sig2_x2),log = TRUE)
  }

  ll0 = sum( apply(dat0,1,f0))
  #ll0=0


  #ll1
  
  f1 <- function(dat){
    dat = as.list(dat)
    counter_bad = 0
    
    constants = log(1/(2*pi*det(sig2_y1y0_x2))) #put back in if not doing top=dnorm()
    integrand <- function(u){
      # use below with constants
      top = diag(exp(-0.5*cbind(u-dat$mu_y0,dat$y1-dat$mu_y1)%*%solve(sig2_y1y0_x2)%*%rbind(u-dat$mu_y0,dat$y1-dat$mu_y1)))
      #top = dmvnorm(cbind(u,dat$y1), mean=cbind(dat$mu_y0,dat$mu_y1), sigma=sig2_y1y0_x2  )
      bmean = dat$mu_y1 + Sig[2,c(1,3),drop=FALSE]%*%solve(Sig_02)%*%rbind(u-dat$mu_y0,dat$x2-dat$mu_x2)
      
      tm2 = proc.time()[1]
      bottom = apply(bmean,2,function(x) pnorm(0,mean=x,sd=sqrt(sig2_y1_y0x2),lower.tail = F))
      #print(paste("bottom time=",proc.time()[1]-tm2))
      
      #ratio = rep(NA,length(top))
      #for(i in 1:length(top)){ratio[i] = as.numeric(top[i]/bottom[[i]][1,1])}
      ratio = top/bottom
      out = as.numeric(ratio)
      
      # if(any(is.nan(out))){browser()}
      if(!all(is.finite(out))){
        counter_bad = 1 + counter_bad
        #print(paste("got a bad one, round ", counter_bad))
        #browser()
        # can't get this to work
        #top = dmvnorm(mpfr(cbind(u,dat$y1),20), mean=cbind(dat$mu_y0,dat$mu_y1), sigma = sig2_y1y0_x2  )
        top = mpfr(log(diag(exp(-0.5*cbind(u-dat$mu_y0,dat$y1-dat$mu_y1)%*%solve(sig2_y1y0_x2)%*%rbind(u-dat$mu_y0,dat$y1-dat$mu_y1)))),100)
        bottom = apply(bmean,2,function(x) pnorm(mpfr(0,100),mean=x,sd=sqrt(sig2_y1_y0x2),lower.tail = F,log.p = TRUE))
        #ratio = rep(NA,length(top))
        #for(j in 1:length(top)){ratio[j] = top[j]/bottom[[j]][1,1]}
        ratio = mapply(function(x,y) as.numeric(exp(x-y)), top, bottom)
        out = unlist(ratio)
      }
      
      
      return(out)
    }
    
    tm1 = proc.time()[1]
    int_result = integrate(integrand,0,Inf)
    #print(paste("int_result=",int_result$value))
    #print(paste("integral time=",proc.time()[1]-tm1))
    
    
    dnorm(dat$x2,mean=dat$mu_x2,sd=sqrt(sig2_x2),log = TRUE) +
      #constants +
      log(int_result$value)
  }
  
  f1_max <- function(dat){
    prec=100
    dat = as.list(dat)
    counter_bad = 0
    
    mean_y1_y0x2 = function(y0,x2, mu_y1, mu_y0, mu_x2){mu_y1 + Sig[2,c(1,3),drop=FALSE]%*%solve(Sig_02)%*%rbind(y0-mu_y0,x2-mu_x2)}
    
    integrand <- function(u){
      # use below with constants
      top_func = function(u){
        mu_y1_y0x2 = mean_y1_y0x2(u,dat$x2, dat$mu_y1, dat$mu_y0, dat$mu_x2)
        dnorm.hack(dat$y1, mu_y1_y0x2, sqrt(sig2_y1_y0x2) ,log=TRUE) + 
          dnorm.hack(u, dat$mu_y0_x2, sqrt(sig2_y0_x2),log=TRUE) + dnorm.hack(dat$x2, dat$mu_x2, sqrt(sig2_x2),log=TRUE)}
      
      bot_func = function(u){
        m = mean_y1_y0x2(u, dat$x2, dat$mu_y1, dat$mu_y0, dat$mu_x2)
        #m = dat$mu_y1 + Sig[2,c(1,3),drop=FALSE]%*%solve(Sig_02)%*%rbind(u-dat$mu_y0,dat$x2-dat$mu_x2)
        pnorm(mpfr(0,prec), m, rep(sqrt(sig2_y1_y0x2),length(m)), lower.tail=FALSE, log.p = TRUE)
      }
      
      #top = new('mpfr', parSapply(clust, u, FUN=top_func))
      top = top_func(u)
      bottom = bot_func(u)
      #top = new("mpfr",sapply(u, FUN=top_func))
      #bottom = new("mpfr", parSapply(clust, u, FUN=bot_func))
      #bottom = new("mpfr", sapply(u, FUN=bot_func))
      #print(cbind(u,top,bottom, top-bottom, exp(top-bottom)))
      out = exp(top-bottom)
      out = as.numeric(out)
      out
    }
    
    tm1 = proc.time()[1]
    int_result = integrate(integrand,0,Inf)
    #print(paste("int_result=",int_result$value))
    #print(paste("integral time=",proc.time()[1]-tm1))
    
    log(int_result$value)
  }
  

  #ll1_old = sum( apply(dat1,1,f1) )
  #ll1 = ll1_old
  ll1 = sum( apply(dat1,1,f1_max) )
  
  print(paste("neg likelihood=",-(ll0+ll1)))
  # assign("JMcounter",JMcounter+1,envir = .GlobalEnv)
  # print(paste("round:",JMcounter))
  print(paste('total time:',proc.time()[1] -tm0))

  #Return the negative log-likelihood, since I will be using optim (which minimizes instead of maximizing)
  -(ll0+ll1)
}

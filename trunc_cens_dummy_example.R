# k: the number of values output can take
# n: the number of samples
# cut: the value at which you censor the outcome

# this file generates a simplified version of the truncation
# you can think of it as a system with flipping coins and rolling dice
# isn't used for anything but checking intuition.

dummy_sim <- function(n,k,cut,N=10000){
  rdu<-function(n,m) sample(1:m,n,replace=T)
  out = matrix(c(rep(NA,length(c(0:(k + 1)))*n)),ncol = length(c(0:(k + 1))))
  n0 = c(rep(NA,n))
  rnd = matrix(c(rep(1,n*N)),ncol = n)
  for(i in 1:n){
    a = rbinom(N,1,.5)
    b = rbinom(N,1,.2+.3*a)
    c = rdu(N,k) + a

    y = ifelse(b==0,0,c)
    n0[i] = length(which((y<cut) & (b>0)))/length(y) #go to second round

    for(j in which((y<cut) & (b>0))){
      m = 2
      while(y[j]<cut){
        ta = rbinom(1,1,.5)
        tb = rbinom(1,1,.2+.3*ta)
        tc = rdu(1,4) + ta

        if((tc>=cut) & (tb>0)){
          a[j] = ta; b[j] = tb; y[j] = tc
          rnd[j,i] = m
        }
        else{m = m+1}
      }
    }
    hist(y)
    out[i,]=c(sapply(c(0:(k+1)),function(x) mean(as.numeric(y==x))))
  }
  list(out,n0,rnd)
}


#P(b=0) = sum_x P(b=0|a=x)P(a=x)
# = (1/2)*(.8) + (1/2)*(1/2)
# = 0.65

#P_gotosecond = P(b=1,c<3,a =x)
# = sum_x P(b=1,c<3|a=x)P(a=x)
# = sum_y sum_x P(b=1,c=y|a=x)P(a=x)
# = P(b=1,c=1|a=0)P(a=0) + P(b=1,c=1|a=1)P(a=1)
# + P(b=1,c=2|a=0)P(a=0) + P(b=1,c=2|a=1)P(a=1)
# by cond'l indep of b&c | a (does not hold for real data)
# = P(b=1|a=0)P(c=1|a=0)P(a=0) + P(b=1|a=1)P(c=1|a=1)P(a=1)
# + P(b=1|a=0)P(c=2|a=0)P(a=0) + P(b=1|a=1)P(c=2|a=1)P(a=1)
# = (1/5)*(1/6)*(1/2) + (1/2)*0*(1/2)
# + (1/5)*(1/6)*(1/2) + (1/2)*(1/6)*(1/2)
# = 0.075

#P_gotothird = P_gotosecond*[P(b=0,c=y,a=x) + P(b=1,c<3,a =x)]
# = P_gotosecond^2 + P_gotosecond*P(b=0,c=y,a=x)
# = 0.005625 + 0.075*sum_x P(b=0|a=x)P(a=x)
# = 0.005625 + 0.075*0.65
# = 0.054375

#P_gotofourth = P_gotothird*[P(b=0,c=y,a=x) + P(b=1,c<3,a=x)]
# = (P_gotosecond^2 + P_gotosecond*P(b=0,c=y,a=x))*[P(b=0,c=y,a=x) + P(b=1,c<3,a=x)]
# = P_gotosecond^3 + P_gotosecond^2*P(b=0,c=y,a=x) + P_gotosecond^2*P(b=0,c=y,a=x) + P_gotosecond*P(b=0,c=y,a=x)


dummy_sim2 <- function(n,k,cut,N=10000){
  require(MASS)
  out = matrix(c(rep(NA,100*n)),ncol = 100)
  rnd = matrix(c(rep(1,n*10000)),ncol = n)
  Sig = matrix(c(1,.2,.3,.2,4,.1,.3,.1,6),ncol = 3,byrow =T)
  for(i in 1:n){
    errs = mvrnorm(N,rep(0,3),Sig)
    a = 3 + errs[,3]
    b = as.numeric(-1 + a + errs[,1]>=0)
    c = 2 + a + errs[,2]

    y = ifelse(b==0,0,c)

    for(j in which((y<cut) & (b>0))){
      m = 2
      while(y[j]<cut){
        e = mvrnorm(1,rep(0,3),Sig)
        ta = 3 + e[3]
        tb = as.numeric(-1 + ta + e[1]>=0)
        tc = 2 + ta + e[2]

        if((tc>=cut) & (tb>0)){
          a[j] = ta; b[j] = tb; y[j] = tc
          rnd[j,i] = m
        }
        else{m = m+1}
      }
    }
    hist(y)
    out[i,]=c(sapply(seq(0,k,length=100),function(x) mean(as.numeric((y>=x) & (y<x+1)))))
  }
  list(out,rnd)
}

dummy_int<-function(u){
  Sigma = matrix(c(1,.4,.4,4),ncol=2,byrow = TRUE)
  mean = matrix(c(-10 + u),ncol = length(u))
  topmean = matrix(cbind(.5 + u, -.3 + u),ncol = length(u),byrow = TRUE)
  temp = apply(mean,2, function(x) pnorm(mpfr(0,100),x,.1,lower.tail = FALSE,log.p = TRUE))
  top = mpfr(log(diag(exp(-0.5*cbind(u-2,rnorm(length(u))-rnorm(length(u)))%*%solve(Sigma)%*%rbind(u-2,rnorm(length(u))-rnorm(length(u)))))),100)
  #top  = apply(rbind(u,topmean),2,function(x) dmvnorm(c(x[1],3),mean=x[2:3],sigma = Sigma,log = TRUE))
  ratio = mapply(function(x,y) as.numeric(exp(x-y)), top, temp)
  return(ratio)
  }

# integrate(dummy_int,0,Inf)


dummy_int2 <-function(u){
  mean = matrix(c(-10 + u),ncol = length(u))
  temp = apply(mean,2, function(x) pnorm(mpfr(0,100),x,.1,lower.tail = FALSE))
  out = mapply(function(x,y) x/y, rep(1,15), temp)
  return(out)
}

dummy_int3 <-function(u){
  mean = matrix(c(-10 + u),ncol = length(u))
  temp = apply(mean,2, function(x) as.bigz(pnorm(0,x,.1,lower.tail = FALSE)))
  out = rep(1,15)/temp
  return(out)
}

testintR<-function(df){

  dummy_int4<-function(u){
    mean = df + u
    top = df - u
    out = top/pnorm(mpfr(0,100),mean = mean,sd=.1,lower.tail = FALSE)
    # ins = as.matrix(cbind(top,mean))
    # out = apply(ins,1,function(x) x[1]/pnorm(mpfr(0,100),x[2],.1,lower.tail = FALSE))
    out
  }

  int_result = integrateR(dummy_int4,0,1e6,ord = 20)
  return(int_result$value)
}


#
# test.data = matrix(c(-10:-8),nrow = 1)
# tosum = apply(test.data,2,testintR)
# Reduce("+",tosum)

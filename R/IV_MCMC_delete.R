likelihood2 <- function(param,type){
  s1 = param[1];s2 = param[2]; t = param[3]; p1 = param[4]; p2 = param[5]; b1 = param[6]; b2 = param[7]

  xpred = p1 + p2*z
  ypred = b1 + b2*xpred
  preds = matrix(c(ypred,xpred),ncol = 2, byrow = F)
  obs = matrix(c(y,x),ncol = 2, byrow = F)
  Sig = matrix(c(s1,t,t,s2),ncol = 2, byrow = T)

  mat = data.frame(y=y,x=x,predy=ypred,predx = xpred)
  singlelikelihoods = apply(mat,1,function(x) dmvnorm(x[1:2],mean = x[3:4],sigma = Sig,log = T))
  sumll1 = sum(singlelikelihoods)

  Sig_err = matrix(c(s1,t,t,s2),ncol = 2, byrow = T)
  A = matrix(c(1,b1,0,1),ncol = 2, byrow = T)
  Sig = A%*%Sig_err%*%t(A)

  mu_y1_x2 = ypred + Sig[1,2]/Sig[2,2]*(x-xpred)
  sig2_y1_x2 = Sig[1,1] - Sig[1,2]^2/Sig[2,2]

  ylik = sum(dnorm(y,mean = mu_y1_x2,sd = sqrt(sig2_y1_x2),log = T))
  xlik = sum(dnorm(x,mean = xpred,sd = sqrt(Sig[2,2]),log = T))
  sumll2 = xlik + ylik

  if(type == 1){return(sumll1)}
  else{return(sumll2)}
}

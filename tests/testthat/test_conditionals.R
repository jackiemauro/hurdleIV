context("test conditional means and variances")
test_that("identity matrices", {
  rho = 0; tau0 = 0; sig_u = 1; tau1 = 0; sig_v = 1
  gamma2 = 0; beta2 = 0;
  x2 = 1.5; y1 = 1.8
  pre = matrix( c(1, rho,    tau0,
                  0,     sig_u^2,    tau1,
                  0,  0, sig_v^2),
                ncol = 3, byrow = T)
  Sig_err = t(pre)%*%pre / ((t(pre)%*%pre)[1,1])
  if(min(eigen(Sig_err)$values)<=0){return(Inf)}
  if((sig_u<=0)|(sig_v<=0)){return(Inf)}

 A = rbind(
    c(1,0,gamma2),
    c(0,1,beta2),
    c(0,0,1)
  )
  Sig = A%*%Sig_err%*%t(A)
  if(min(eigen(Sig)$values)<=0){return(Inf)}

  mu_y0 = 1
  mu_y1 = 2
  mu_x2 = 3

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
  sig2_y1y0_x2 = Sig[1:2,1:2] - Sig[1:2,3]%*%solve(sig2_x2)%*%Sig[3,1:2]
  sig2_y1y0_x2[upper.tri(sig2_y1y0_x2)] <- sig2_y1y0_x2[lower.tri(sig2_y1y0_x2)]
  if(any(eigen(sig2_y1y0_x2)$value<0)){return(-Inf)}

  #Parameters for y1star and y0star given x2
  Sig_02 = Sig[-2,-2]
  sig2_y1_y0x2 = Sig[2,2] - Sig[2,c(1,3),drop=FALSE]%*%solve(Sig_02)%*%Sig[c(1,3),2,drop=FALSE]

  allsame <- function(x,tol=0){abs(max(x)-min(x))<=tol}
  expect_true(allsame(c(mu_y0,mu_y0_x2,mu_y0_y1x2)))
  expect_true(allsame(c(mu_y1,mu_y1_x2)))
  expect_true(allsame(c(sig2_y0_x2,sig2_y0_y1x2,sig2_y0_y1x2,Sig[1,1])))
  expect_true(allsame(c(sig2_y1_x2,sig2_y1_y0x2,Sig[2,2])))
  expect_true(all(sig2_y1y0_x2 == Sig[1:2,1:2]))

})


context("test conditional means and variances")
test_that("easy 3x3", {
  rm(list = ls())
  sig_u = 2; sig_v = 3
  x2 = 1.5; y1 = 1.8
  Sig = matrix( c(1, .2,    2,
                  .2,     sig_u^2,    .1,
                  2,  .1, sig_v^2),
                ncol = 3, byrow = T)

  if(min(eigen(Sig)$values)<=0){return(Inf)}

  mu_y0 = 1
  mu_y1 = 2
  mu_x2 = 3

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
  sig2_y1y0_x2 = Sig[1:2,1:2] - Sig[1:2,3]%*%solve(sig2_x2)%*%t(Sig[1:2,3])
  sig2_y1y0_x2[upper.tri(sig2_y1y0_x2)] <- sig2_y1y0_x2[lower.tri(sig2_y1y0_x2)]
  if(any(eigen(sig2_y1y0_x2)$value<0)){return(-Inf)}

  #Parameters for y1star and y0star given x2
  Sig_02 = Sig[-2,-2]
  sig2_y1_y0x2 = Sig[2,2] - Sig[2,c(1,3),drop=FALSE]%*%solve(Sig_02)%*%Sig[c(1,3),2,drop=FALSE]

  allsame <- function(x,tol=0){abs(max(x)-min(x))<=tol}
  expect_equal(mu_y0_x2, 2/3); expect_equal(sig2_y0_x2,1-(4/9))
  expect_equal(mu_y1_x2,1.983333,tol=1e-6); expect_equal(sig2_y1_x2,3.998889,tol=1e-6)
  expect_equal(sig2_y1_y0x2[1,1],3.942)
  expect_equal(sig2_y1y0_x2,matrix(c(1,.2,.2,4),ncol=2,byrow=T)-(1/9)*matrix(c(4,.2,.2,.01),ncol=2,byrow=T))
  expect_equal(mu_y0_y1x2,0.6585163,tol = 1e-6)
  expect_equal(sig2_y0_y1x2[1,1],.5476521,tol = 1e-6)

})

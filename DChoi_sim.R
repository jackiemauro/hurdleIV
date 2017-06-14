#### this is to test DChoi's likelihood (or my interpretation of it) ###
# source("loglik_cragg_dChoi.R")
# source("cragg_errors_dChoi.R")
require(mvtnorm)

# exogenous variable/instrument
N = 10000
x1 = rnorm(N,3,1)
z = rnorm(N,.2,1)

# parameters
cov = matrix(c(1,.2,.1,.2,4,.3,.1,.3,6),ncol = 3,byrow = T)
pi = c(1,2,1)
gamma = c(.1,.2,-.1)
beta = c(.3,.1,-.2)

# get errors and generated variables
dat = cragg_errs_DC(cov=cov,pi=pi,gamma=gamma,beta=beta,x1=x1,z=z,n=N)
y1 = c(dat$yStar)
x2 = c(dat$endog)

# start values and true values
start = c(6,4,.1,.1,.1,.1,.1,.1,.1,.1,.1,1,1,1)
chol = chol(cov)
true = c(chol[3,3],chol[2,2],chol[2,3],chol[1,3],chol[1,2],beta,gamma,pi)

# run the optimizer and examine
ll = optim(start,loglik_dChoi,hessian = T)
ll = optim(true,loglik_dChoi,hessian = T)
nm = c('sigv','sigu','t0','t1','rho','beta11','beta12','beta2','gamma11','gamma12','gamma2','pi11','pi12','pi2')
data.frame(name = nm,est = ll$par,true = true) #how are the parameters
min(eigen(ll$hessian)$values); max(eigen(ll$hessian)$values) #should both be >0

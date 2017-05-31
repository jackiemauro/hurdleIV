#### this is to test MG's likelihood (or my interpretation of it) ###

# exogenous variable/instrument
N = 1000
x1 = rnorm(N,3,1)
z = rnorm(N,.2,1)

# parameters
cov = matrix(c(1,.2,.1,.2,4,.3,.1,.3,6),ncol = 3,byrow = T)
pi = c(1,2,1)
gamma = c(.1,.2,-.1)
beta = c(.3,.1,-.2)

# get errors and generated variables
dat = cragg_errs3(cov=cov,pi=pi,gamma=gamma,beta=beta,x1=as.matrix(x1),z=as.matrix(z),n=N)
y1 = c(dat$yStar)
x2 = c(dat$endog)

# start values and true values
start = c(6,4,.1,.1,.1,.1,.1,.1,.1,.1,.1,1,1,1)
chol = chol(cov)
true = c(chol[3,3],chol[2,2],chol[2,3],chol[1,3],chol[1,2],beta,gamma,pi)

# run the optimizer and examine
#ll = optim(start,loglik_dChoi,hessian = T)
ll = optim(true,loglik_MG,hessian = T)
cbind(ll$par,true) #how are the parameters
min(eigen(ll$hessian)$values); max(eigen(ll$hessian)$values) #should both be >0

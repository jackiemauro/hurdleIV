#### this is to test MG's likelihood (or my interpretation of it) ###

# exogenous variable/instrument
N = 1000
x1 = rnorm(N,3,1)
z = rnorm(N,.2,1)

# parameters
cov = matrix(c(1,.2,.1,.2,4,.3,.1,.3,6),ncol = 3,byrow = T)
pi = c(1,2,1)
gamma = c(.1,.2,-.1)
beta = c(.3,-1,-.2)

# get errors and generated variables
dat = cragg_errs_MG(cov=cov,pi=pi,gamma=gamma,beta=beta,x1=as.matrix(x1),z=as.matrix(z),n=N)
y1 = c(dat$yStar)
x2 = c(dat$endog)

# what if i sample split?
sample = sample(c(1:N),N/2)
df1 = data.frame(y1 = y1, x2 = x2, z = z, x1 = x1)[sample,]
df2 = data.frame(y1 = y1, x2 = x2, z = z, x1 = x1)[-sample,]

x2reg = lm(x2~x1+z,data = df1)
x2hat = predict(x2hat,newdata = df2)
probit = glm(as.numeric(y1>0)~x2hat + x1, data = df2, family = binomial(link = probit))
linear = lm(y1~x1 + x2hat, data = df2)

# start values and true values
start = c(6,4,.1,.1,.1,.1,.1,.1,.1,.1,.1,1,1,1)
chol = chol(cov)
true = c(chol[3,3],chol[2,2],chol[2,3],chol[1,3],chol[1,2],beta,gamma,pi)

# run the optimizer and examine
#ll = optim(start,loglik_MG,hessian = T)
JMcounter = 0
ll = optim(true,loglik_MG,hessian = T)
nms = c("sigv", "sigu","tau1","tau0","rho","beta11","beta12","beta2","gamma11","gamma12","gamma2","pi11","pi12","pi2")
View(cbind(nms,ll$par,true)) #how are the parameters
min(eigen(ll$hessian)$values); max(eigen(ll$hessian)$values) #should both be >0

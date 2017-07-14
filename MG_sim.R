# #### this is to test MG's likelihood (or my interpretation of it) ###
#
# # exogenous variable/instrument
# N = 1000
# x1 = rnorm(N,3,1)
# z = rnorm(N,.2,1)
#
# parameters
cov = matrix(c(1,.2,.1,.2,4,.3,.1,.3,6),ncol = 3,byrow = T)
pi = c(1,2,1)
gamma = c(.1,.2,-.1)
beta = c(.3,-1,-.2)
#
# # get errors and generated variables
# dat = cragg_errs_MG(cov=cov,pi=pi,gamma=gamma,beta=beta,x1=as.matrix(x1),z=as.matrix(z),n=N)
# y1 = c(dat$yStar)
# x2 = c(dat$endog)
# write.csv(file = 'testdata.csv',data.frame(y1,x1,x2,z))

# # what if i sample split?
# sample = sample(c(1:N),N/2)
# df1 = data.frame(y1 = y1, x2 = x2, z = z, x1 = x1)[sample,]
# df2 = data.frame(y1 = y1, x2 = x2, z = z, x1 = x1)[-sample,]
#
# x2reg = lm(x2~x1+z,data = df1)
# x2hat = predict(x2reg,newdata = df2)
# probit = glm(as.numeric(y1>0)~x2hat + x1, data = df2, family = binomial(link = probit))
# linear = lm(y1~x1 + x2hat, data = df2)

dat <- read.csv('testdata.csv')[,-1]
y1 <- dat$y1; x1<- dat$x1; x2<-dat$x2; z <- dat$z

# start values and true values
start = c(6,4,.1,.1,.1,.1,.1,.1,.1,.1,.1,1,1,1)
chol = chol(cov)
true = c(chol[3,3],chol[2,2],chol[2,3],chol[1,3],chol[1,2],beta,gamma,pi)

# run the optimizer and examine
JMcounter = 0
ll = optim(true,loglik_MG,hessian = T)
# write.csv(file = "hessian.csv", ll$hessian)
# write.csv(file = "parametersOut.csv",ll$par)

# after making sig2_y0_y1x2 correction
write.csv(file = "hessianNew.csv", ll$hessian)
write.csv(file = "parametersOutNew.csv",ll$par)

#ll = optim(true,loglik_MG,hessian = T)
nms = c("sigv", "sigu","tau1","tau0","rho","beta11","beta12","beta2","gamma11","gamma12","gamma2","pi11","pi12","pi2")
View(cbind(nms,ll$par,true)) #how are the parameters
min(eigen(ll$hessian)$values); max(eigen(ll$hessian)$values) #should both be >0

bad_eig <- eigen(ll$hessian)$vec[,14]

#####check along bad eigens ####
hess <- as.matrix(read.csv('hessian.csv')[,-1])
pars <- as.matrix(read.csv('parametersOut.csv')[,-1])
eigs <- eigen(hess)
bad_eig <- matrix(eigs$vectors[,14],ncol = 1)
dat <- read.csv('testdata.csv')[,-1]
attach(dat)

#create dataframe the varies pars along bad eig's dimension
n = 3; seq = seq(-n,n, by =.25 )
JMcounter = 0
out=sapply(seq,function(x) loglik_MG(pars+bad_eig*(x/n)))
plot(seq/n,out,xlab = "deviation",ylab="objective",pch = 19)

# check first and second differences
first = diff(out); fpos = which(first>0); first_pos = first[fpos]; first_neg = first[-fpos]
sec = diff(out,2); spos = which(sec>0); sec_pos = sec[spos]; sec_neg = sec[-spos]

par(mfrow = c(1,2))
seq1 = seq[-1]; seq2 = seq[-c(1,2)]
plot(seq1/n,diff(out),xlab = "deviation",ylab="difference",main = "first difference",type = 'n')
points((seq1/n)[fpos],first_pos,col = 'black',pch = 19);points(seq1[-fpos]/n,first_neg,col = 'red',pch = 19)
plot(seq2/n,diff(out,2),xlab = "deviation",ylab="difference",main = "second difference",type = 'n')
points((seq2/n)[spos],sec_pos,col = 'black',pch = 19);points(seq2[-spos]/n,sec_neg,col = 'red',pch = 19)
par(mfrow = c(1,1))

# see how close you can get with regression
x = seq(1:14); y = eigs$values
plot(x, y, pch = 19, ylab = "eigenvalues")
text(x, y, round(y, 2), cex=0.5)

reg1 <- lm( (true - pars) ~ -1+eigs$vectors[,-c(1:7)] )
reg2 <- lm( (true - pars) ~ -1+eigs$vectors[,c(13:14)] )

fit1 <- pars + fitted(reg1)
fit2 <- pars + fitted(reg2)

compare = data.frame(true = true, estimated = pars, last2 = fit2, last7 = fit1)
write.csv(compare,file = "getThetafromReg.csv")

# compare r squared across regressions
f <- function(x){summary(lm( (pars - true) ~ eigs$vectors[,-c(1:x)] ))$r.squared}
rs <- sapply(c(1:13),f)
plot(c(1:13), rs, xlab = 'eigenvalue rank', ylab = 'r squared',xlim=c(0,15),pch=19)
text(c(1.5:13.5), rs-.01, round(eigs$values[2:14]), cex=0.6)

# compare residual standard error across regressions
f <- function(x){summary(lm( (pars - true) ~ eigs$vectors[,-c(1:x)] ))$sigma}
sig <- sapply(c(1:13),f)
plot(c(1:13), sig, xlab = 'eigenvalue rank', ylab = 'residual SE',xlim=c(0,15),pch=19)
text(c(1.5:13.5), sig, round(eigs$values[2:14]), cex=0.6)

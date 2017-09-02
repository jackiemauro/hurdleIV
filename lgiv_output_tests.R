

# does really well, except on tau's
dat = hurdle.IV.sim()
out = hurdle.IV(y~exog1 + endog,
                inst = inst1,
                endog = endog,
                exog = exog1,
                data = dat)

# small sample
dat = hurdle.IV.sim(n=50)
out = hurdle.IV(y~exog1 + endog,
                inst = inst1,
                endog = endog,
                exog = exog1,
                data = dat)

# same as above
dat = hurdle.IV.sim(beta = c(1,2,3))
out = hurdle.IV(y~exog1 + endog,
                inst = inst1,
                endog = endog,
                exog = exog1,
                data = dat)

# actually maybe a little better?
dat = hurdle.IV.sim(tau0=-1,tau1 = 3, rho = 1)
out = hurdle.IV(y~exog1 + endog,
                inst = inst1,
                endog = endog,
                exog = exog1,
                data = dat)

#multiple covariates
bs = c(.05,.06,-.01,.01,.02)
gs = c(-.2,.08,.3,-.3,.07)
ps = c(1,-1,2,2,3)
mn = c(1,1,-1)
s = c(1,1,1)

dat = hurdle.IV.sim(beta=bs, gamma = gs, pi = ps,
                    exog_mean = mn, exog_sd = s)
out = hurdle.IV(y~exog1 + exog2 + exog3 + endog,
                inst = inst1,
                endog = endog,
                exog = c(exog1,exog2,exog3),
                data = dat)

#multiple covariates
bs = c(.05,.06,-.01,.01,.02)
gs = c(-.2,.08,.3,-.3,.07)
ps = c(1,-1,2,2,3)
mn = c(1,1,-1)
s = c(1,1,1)

dat = hurdle.IV.sim(beta=bs, gamma = gs, pi = ps,
                    exog_mean = mn, exog_sd = s)
out = hurdle.IV(y~exog1 + exog2 + exog3 + endog,
                inst = inst1,
                endog = endog,
                exog = c(exog1,exog2,exog3),
                data = dat)

#multiple covariates, one binary (just to check if runs)
bs = c(.05,.06,-.01,.01,.02)
gs = c(-.2,.08,.3,-.3,.07)
ps = c(1,-1,2,2,3)
mn = c(1,1,-1)
s = c(1,1,1)

dat = hurdle.IV.sim(beta=bs, gamma = gs, pi = ps,
                    exog_mean = mn, exog_sd = s)
dat = as.data.frame(cbind(dat,exog4 = rbinom(length(dat[,1]),1,.3)))
out = hurdle.IV(y~exog1 + exog2 + exog3 + exog4 + endog,
                inst = inst1,
                endog = endog,
                exog = c(exog1,exog2,exog3,exog4),
                data = dat)

#multiple covariates
bs = c(.05,.06,-.01,.01,.02)
gs = c(-.2,.08,.3,-.3,.07)
ps = c(1,-1,2,2,3)
mn = c(1,1,-1)
s = c(1,1,1)

dat = hurdle.IV.sim(beta=bs, gamma = gs, pi = ps,
                    exog_mean = mn, exog_sd = s)
out = hurdle.IV(y~exog1 + exog2 + exog3 + endog,
                inst = inst1,
                endog = endog,
                exog = c(exog1,exog2,exog3),
                data = dat)

#multiple covariates, missing data
bs = c(.05,.06,-.01,.01,.02)
gs = c(-.2,.08,.3,-.3,.07)
ps = c(1,-1,2,2,3)
mn = c(1,1,-1)
s = c(1,1,1)

dat = hurdle.IV.sim(beta=bs, gamma = gs, pi = ps,
                    exog_mean = mn, exog_sd = s)
samp = sample(1:70000, 300, replace = FALSE); matdat = as.matrix(dat)
matdat[samp]<-NA
dat = as.data.frame(matdat)

out = hurdle.IV(y~exog1 + exog2 + exog3 + endog,
                inst = inst1,
                endog = endog,
                exog = c(exog1,exog2,exog3),
                data = dat)

#multiple covariates--one factor (just to see if it runs)
bs = c(.05,.06,-.01,.01,.02)
gs = c(-.2,.08,.3,-.3,.07)
ps = c(1,-1,2,2,3)
mn = c(1,1,-1)
s = c(1,1,1)

dat = hurdle.IV.sim(beta=bs, gamma = gs, pi = ps,
                    exog_mean = mn, exog_sd = s)
dat = as.data.frame(cbind(dat,exog4 = sample(1:500,dim(dat)[1],replace = T)))
out = hurdle.IV(y~exog1 + exog2 + exog3 + as.factor(exog4) + endog,
                inst = inst1,
                endog = endog,
                exog = c(exog1,exog2,exog3,exog4),
                data = dat)


#MCMC
bs  = c(0.05, 0.06, 0.02); gs = c(-0.2, 0.8, 0.07); pis = c(1, -1, 3)
sdy1 = 2; sdx2 = 5; r = .2; t0 = .3; t1 = .1
dat = hurdle.IV.sim()
out = hurdle.IV.MCMC(y~exog1 + endog,
                inst = inst1,
                endog = endog,
                exog = exog1,
                data = dat,
                options = list(cholesky = T
                               , maxit = 5000
                               , trace = 0
                               , method = "BFGS"),
                k = 1000)

nms = c('rho','sig2y1','t0','t1','sig2x2','beta1','beta2','beta3','gamma1','gamma2','gamma3','pi1','pi2','pi3')
true = c(r,sdy1,t0,t1,sdx2,bs,gs,pis)
mns = apply(out,2,mean);sds = apply(out,2,sd)
View(cbind(true,mns))
for(i in 1:14){hist(out[,i],main = nms[i]);abline(v = true[i],col = 'red');abline(v = mns[i],col = 'green')}

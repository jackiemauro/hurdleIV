

# does really well, except on tau's
dat = hurdle.IV.sim()
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

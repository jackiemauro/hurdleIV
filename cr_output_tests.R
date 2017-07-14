devtools::load_all();detach(pars);detach(mf)

dat = hurdle.IV.sim(tau0=-1,tau1 = 3, rho = 1,type='cragg')
out = hurdle.IV(y~exog1 + endog,
                inst = inst1,
                endog = endog,
                exog = exog1,
                type = 'cragg',
                data = dat)

dat = hurdle.IV.sim(pi = c(1,-1,2,3)
                    ,gamma = c(-.2,.8,.1,.07)
                    ,beta = c(.1,.01,.02,.02)
                    ,exog_mean = c(1,2)
                    ,exog_sd = c(1,1)
                    ,type='cragg')
out = hurdle.IV(y~exog1 + exog2 + endog,
                inst = inst1,
                endog = endog,
                exog = c(exog1,exog2),
                type = 'cragg',
                data = dat)

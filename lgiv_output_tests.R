

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

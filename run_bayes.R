library('rjags')

jags <- jags.model('model_priors.bug',
                   data = list('x11' = rep(1,length(dat$exog1)),
                               'x12' = dat$exog1,
                               'x21' = dat$endog,
                               'z1' = dat$inst1,
                               'y' = dat$y,
                               'N' = dim(dat)[1]),
                   n.chains = 4,
                   n.adapt = 100)

update(jags, 1000)

jags.samples(jags,
             c('a', 'b'),
             1000)

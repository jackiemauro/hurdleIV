library('rjags')

# see: https://github.com/johnmyleswhite/JAGSExamples/commit/ad1e8716607d50f0eee62411bb69183dd93e5f99

# simulate data
N=5000
x = rnorm(N,3,1)
a = -2; b = 2
y = a*x + b + rnorm(N,0,2)
y = y*as.numeric(y>=0)
dataList = list( y = y , x = x, N = N )

# model
jags <- jags.model('example.bug',
                   data = dataList,
                   n.chains = 4,
                   n.adapt = 100)

update(jags, 1000)

jags.samples(jags,
             c('a', 'b','sigma'),
             1000)




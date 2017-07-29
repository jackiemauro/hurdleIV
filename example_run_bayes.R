library('rjags')

# see: http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/Programs/JagsCensoringExample.R

# THE MODEL
modelstring = "
model
{
  for (i in 1:N)
  {
  isCensored[i] ~ dinterval( y[i] , censorLimitVec[i] )
  y[i] ~ dnorm(mu[i], tau)
  mu[i] <- a * x[i] + b
  }

  a ~ dnorm(0, 0.0001)
  b ~ dnorm(0, 0.0001)

  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 10000)
  }

"
writeLines(modelstring,con="model.txt")

# simulate data
N=1000
x = rnorm(N,3,1)
a = -2; b = 1
y = a + b*x + rnorm(N,0,2)
censorLimitVec = rep( 0 , length(y) )
isCensored = ( y <= censorLimitVec )
y[isCensored] <- NA
dataList = list( y = y , x = x, N = N
                 , isCensored = as.numeric(isCensored)
                 , censorLimitVec = censorLimitVec
)

# INTIALIZE THE CHAINS.
sigmaInit = 15
aInit = 100; bInit = 100
# intial values of censored data:
yInit=y
yInit[]=NA
yInit[isCensored] = censorLimitVec[isCensored]+1
initsList = list(sigma = sigmaInit, a = aInit, b = bInit, y=yInit )

# model
parameters = c( "a" , "b" , "sigma" )
adaptSteps = 100 # overkill, just shows general form
burnInSteps = 2000 # overkill, just shows general form
nChains = 3
numSavedSteps=50000
thinSteps=1
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

jags <- jags.model('model.txt',
                   data = dataList,
                   inits = initsList,
                   n.chains = nChains,
                   n.adapt = adaptSteps)
codaSamples = coda.samples( jags , variable.names=parameters ,
                            n.iter=nIter
                            , thin=thinSteps
                            )

checkConvergence = FALSE
if ( checkConvergence ) {
  openGraph(width=7,height=7)
  autocorr.plot( codaSamples[[1]] , ask=FALSE )
  show( gelman.diag( codaSamples ) )
  effectiveChainLength = effectiveSize( codaSamples )
  show( effectiveChainLength )
}


# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )

par( mar=c(3.5,1.5,3.0,0.5) , mgp=c(2,0.7,0) )
desiredMean = a*3+b; desiredSD = 2
# plot the data:
xLim = desiredMean+c(-3.5*desiredSD,3.5*desiredSD)
xBreaks = seq( xLim[1] , xLim[2] , desiredSD/2 )
xGrid = seq( xLim[1] , xLim[2] , length=(N+1) )
histInfo = hist( y , main="Data" , xlab="y" , freq=F ,
                 xlim=xLim , breaks=c(min(y,na.rm=T),xBreaks) ,
                 col="grey" , border="white" )
lines( xGrid , dnorm( xGrid , mean=desiredMean , sd=desiredSD )/pnorm(1) , lwd=3 )
abline(v=censorLimit,lty="dotted",col="red")
# plot the posterior:
histInfo = plotPost( mcmcChain[,"a"] , main="Mu" , xlab=expression(mu) ,
                     showMode=FALSE )
histInfo = plotPost( mcmcChain[,"sigma"] , main="Sigma" , xlab=expression(sigma) ,
                     showMode=TRUE )

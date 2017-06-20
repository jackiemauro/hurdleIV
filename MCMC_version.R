# from : https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/

prior <- function(params){
  sig_v = params[1]
  sig_u = params[2]
  tau1 = params[3]
  tau0 = params[4]
  rho = params[5]
  beta11 = params[6]
  beta12 = params[7]
  beta2 = params[8]
  gamma11 = params[9]
  gamma12 = params[10]
  gamma2 = params[11]
  pi11 = params[12]
  pi12 = params[13]
  pi2 = params[14]

  sigvp = dunif(sig_v, min=0, max=100, log = T)
  sigup = dunif(sig_u, min=0, max=100, log = T)
  tau1p = dunif(tau1, min = 0, max = 1)
  tau0p = dunif(tau0, min = 0, max = 1)
  rhop = dunif(rho, min = 0, max = 1)
  coefp = sapply(c(beta11,beta12,beta2,gamma11,gamma12,gamma2,pi11,pi12,pi2),
                function(x) dnorm(x,mean=0,sd=10))
  return(sigvp+sigup+tau1p+tau0p+rhop+sum(coefp))
}

posterior <- function(param){
  ll = loglik_MG(param); pr = prior(param)
  if(is.nan(ll)){cat('ll=',ll)}
  if(is.infinite(ll)){print(paste('post returning Inf')); return(Inf)}
  else{print(paste('post returning',ll+pr));return (ll + pr)}
}

proposalfunction <- function(param){
  return(rnorm(14,mean = param, sd= c(rep(.1,14))))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,14))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])

    post1 = posterior(proposal); post2 = posterior(chain[i,])
    print(paste("post1=",post1))
    print(paste("post2=",post2))
    probab = exp(post1 - post2)
    if(is.infinite(post1) & is.infinite(post2)){probab = 0}
    print(paste("probab=",probab))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

startvalue = c(rep(.1,14))
chain = run_metropolis_MCMC(startvalue, 10000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

myMCMC<-function(params,k){
  likelihood<-loglik_lgiv

  prior <- function(params){
    #need to get this to be wishart
    invsds = dnorm(params[1:len_cov],sd = 100,log = T)
    sds = 1/invsds
    lasts = Reduce("+",lapply(params[(len_cov+1):length(params)],function(x) dnorm(x,100,log=T)))
    lasts + sum(sds)
  }

  posterior <- function(param){
    return (likelihood(param) + prior(param))
  }

  proposalfunction <- function(param){
    return(rnorm(length(param),mean = param, sd= rep(1,length(param))))
  }

  run_my_MCMC <- function(startvalue, iterations){
    chain = array(dim = c(iterations+1,length(startvalue)))
    chain[1,] = startvalue
    for (i in 1:iterations){
      proposal = proposalfunction(chain[i,])

      probab = exp(posterior(proposal) - posterior(chain[i,]))
      if ((!is.na(probab)) & runif(1) < probab){
        chain[i+1,] = proposal
      }else{
        chain[i+1,] = chain[i,]
      }
    }
    return(chain)
  }

  out = run_my_MCMC(params, k)
  return(out)
}


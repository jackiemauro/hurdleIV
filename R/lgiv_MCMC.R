myMCMC<-function(params,dat,k,pars,mcmcAlgo = "HARM", thin = 1,...){
  model <- function(parm, data) {
    log_lik<-loglik_lgiv(parm)

    pieces = name.pieces(parm)

    covs = pieces$cov_start
    sigy1 = exp(covs[2])
    sigx = exp(tail(covs,1))
    log_prior <- sum(dnorm(c(pieces$pi[[1]], pieces$gamma, pieces$beta, covs),0,100,log = T))+dlnorm(sigy1, 0, 4, log = T)+dlnorm(sigx, 0, 4, log = T)

    post_lik <- log_lik + log_prior

    # This list is returned and has to follow a the format specified by
    # LaplacesDemon.
    list(LP = post_lik, Dev = -2 * log_lik, Monitor = c(post_lik, sigx,sigy1), yhat = NA,parm = parm)
  }

  pieces = name.pieces(params)
  # betnames = sapply(c(1:length(pieces$beta)),function(x) paste('beta',x,sep = ""))
  # gamnames = sapply(c(1:length(pieces$gamma)),function(x) paste('gamma',x,sep = ""))
  # pinames = sapply(c(1:length(pieces$pi[[1]])),function(x) paste('pi',x,sep = ""))
  names = c('rho','logsigy1','tau0','tau1','logsigx',names(pieces$beta),names(pieces$gamma),names(pieces$pi[[1]]))

  # Running the model
  data_list <- list(N = length(dat$outcome), y = dat$outcome,x = dat$endog, z = dat$inst1,
                    mon.names = c("post_lik", "sigx","sigy1"),
                    parm.names = names)
  mcmc_samples <- LaplacesDemon(Model = model, Data = data_list, Iterations = k,
                                Algorithm = mcmcAlgo, Thinning = thin,Initial.Values = params,...)

  return(mcmc_samples)
}








myMCMC2<-function(params,dat,k,pars){
  #uses own MCMC
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


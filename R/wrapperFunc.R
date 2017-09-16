#' Hurdle IV regression
#'
#' @description Run a hurdle IV regression, either lognormal IV or cragg IV.
#' You should specify the model you want to fit. Currently this does not have
#' a nice summary function for the output.
#'
#' @param formula the second stage regression: y~exogenous + endogenous
#' @param inst a vector of your instrument(s): c(inst1,inst2)
#' @param endog a vector of your endogenous variable(s): c(end1,end2)
#' @param exog a vector of your exogenous variable(s): c(ex1,ex2)
#' @param data the dataframe
#' @param endog_reg a list of endogenous regression formulae. By default will
#' estimate endogenous ~ all exogenous variables + all instruments
#' @param start_val an optional list of start values. By default the function
#' will find start values using simple linear/probit regressions of the specified
#' formulae
#' @param type either "lognormal" or "cragg"
#' @param options: cholesky true or false; maxit = maximal number of iterations;
#' trace = 0; method = "BFGS". Options are similar to those for optim, see
#' optim documentation for more details
#'
#' @return Returns estimated parameters as well as a the hessian and standard
#' deviations


hurdle.IV<-function(formula,
                    inst,
                    endog,
                    exog,
                    data,
                    endog_reg = list(),
                    start_val = list(),
                    type = "lognormal",
                    options = list(cholesky = T
                                   , maxit = 5000
                                   , trace = 0
                                   , method = "BFGS")
                    ){

  # helper fn
  rename.input <- function(input){
    if(class(input)=="name"){
      name = toString(input)
    }
    else if(class(input) == "call"){
      text = input[-1]
      k = length(text)
      name = NULL
      for(i in 1:k){name[i] = toString(text[[i]])}
      return(name)
    }
  }

  ###### format your data, grouping variable types #####
  data = data[complete.cases(data),]
  attach(data)
  inst_mat = as.data.frame(matrix(inst,nrow = dim(data)[1], byrow = F))
  inst_names = rename.input(substitute(inst))
  names(inst_mat)<-inst_names

  endog_mat = as.data.frame(matrix(endog, nrow = dim(data)[1], byrow = F))
  endog_names = rename.input(substitute(endog))
  names(endog_mat)<- endog_names

  if(is.null(exog)){
    # if you want to run with no covariates
    exog_mat = data.frame(none = rep(0,dim(data)[1]))
    exog_names = NULL
  }
  else{
    exog_mat = as.data.frame(matrix(exog, nrow = dim(data)[1], byrow = F))
    exog_names = rename.input(substitute(exog))
  }
  names(exog_mat)<- exog_names


  outcome = eval(parse(text=paste("data$",formula[[2]],sep = "")))
  y_mat = model.matrix(formula)
  y_colnames = colnames(y_mat)
  y_exog = exog_mat[names(exog_mat) %in% y_colnames]
  y_endog = endog_mat[names(endog_mat) %in% y_colnames]
  if(any(names(inst_mat)%in% y_colnames)){
    print("Error: main regression cannot include instruments")
    return(NA)
  }
  ######### checks #########
  #check you have a formula
  tryCatch({
    form = as.formula(formula)
  }, error = function(e){
    detach(data)
    stop("Error: formula must be a formula, eg: y~x1+x2")
  }
  )
  # check the formula includes endogenous variables
  if(dim(y_endog)[2]==0){
    detach(data)
    stop("Error: endogenous variables must be included in main regression")
  }

  # check no more endogenous variables than instruments
  if(length(inst)<length(endog)){
    detach(data)
    stop("Error: More endogenous variables than instruments")
  }

  # if not specified, replace endog_reg with formula that has all exog's, all inst's
  if(!is.list(endog_reg)){
    stop("Error: endog_reg must be a list (can be empty)")
  }




  ##### make endog reg #####
  if(length(endog_reg) == 0){
    a = names(endog_mat)[1]
    d = paste(names(inst_mat), collapse = "+")
    if(is.null(exog)){endog_reg = list(as.formula(paste(a,"~",d)))}
    else{
      b = paste(names(exog_mat),collapse = "+")
      endog_reg = list(as.formula(paste(a,"~",b,"+",d)))
    }
  }
  else{
    if(length(endog_reg)!=dim(endog_mat)[2]){
      detach(data)
      stop("Error: Need one regression for each endogenous variable")
    }
  }

  ER_mat = lapply(endog_reg, function(x) model.matrix(x))
  ER_colnames = lapply(ER_mat, function(x) colnames(x))
  ER_exog = lapply(ER_colnames, function(x) exog_mat[names(exog_mat) %in% x])
  ER_inst = lapply(ER_colnames, function(x) inst_mat[names(inst_mat) %in% x])

  ### drop NA's
  if(is.null(exog)){mf = cbind(outcome, inst_mat, endog_mat)}
  else{mf = cbind(outcome, exog_mat, inst_mat, endog_mat)}
  mf = mf[complete.cases(mf),] #drop NA's
  detach(data)
  attach(mf)

  ############# get start values #######
  # if start values aren't specified, get start values
  if(length(start_val) == 0){
    start_val = start.val(formula = update(form,outcome~.)
                          , endog_reg = endog_reg
                          , data = mf
                          , type = type)
  }

  else{
    if(class(start_val) != "list"){
      detach(mf)
      stop("Error: start values must be a list")
    }
    if(length(start_val[['beta']])!=length(start_val[['gamma']])){
      detach(mf)
      stop("Error: beta and gamma must be equal length vectors")
    }
    k = dim(endog_mat)[2]
    if(length(start_val$pi)!=k | length(start_val$tau0)!= k | length(start_val$tau1)!=k){
      detach(mf)
      stop("Error: tau0, tau1 and number of pi coordinates must equal number of endogenous variables")
    }
    if(is.null(names(start_val$gammas))){
      names(start_val$gamma)<-c("Intercept", exog_names, endog_names)
    }
    if(is.null(names(start_val$beta))){
      names(start_val$beta)<-c("Intercept", exog_names, endog_names)
    }
  }

  # start cov values should be cholesky-fied if that option is true
  cov_start = make.cov(rho=start_val$rho
                       ,tau0=start_val$tau0
                       ,tau1=start_val$tau1
                       ,y_sd=start_val$y_sd
                       ,endog_sd=start_val$endog_sd)
  if(options$cholesky == T){
    cov_start = chol(cov_start)
  }


  ########## run optimizer #############
  # save info about parameters and model matrices
  pars = list(len_cov = length(which(upper.tri(cov_start, diag = T)))-1
              ,num_endog = dim(endog_mat)[2]
              ,num_betas = dim(y_mat)[2]
            #  ,num_pis = lapply(ER_mat, function(x) dim(x)[2])
              ,numpis = lapply(ER_mat, function(x) dim(x)[2])
              ,myChol = options$cholesky
              ,ER_mat = ER_mat
              ,y_mat = y_mat
              ,endog_mat = endog_mat
              ,inst_names = inst_names
              ,exog_names = exog_names
              ,endog_names = endog_names)
  attach(pars)

  cov_in = cov_start[upper.tri(cov_start, diag = T)][-1]
  start_vec = c(cov_in, start_val$beta, start_val$gamma, unlist(start_val$pi))

  if(type == "lognormal"){
    ll = loglik_lgiv
  }
  else if(type == "cragg"){
    ll = loglik_craggiv
  }

  out = optim(start_vec
              , ll
              , method= options$method
              , hessian  = T
              , control = list(maxit = options$maxit,trace = options$trace)
  )
  name_pars = name.pieces(out$par)
  name_pars$cov = reconstitute.cov(vals=name_pars$cov_start
                                   ,num=num_endog
                                   ,chol=myChol)
  name_pars$cov_start <-NULL
  detach(mf)
  detach(pars)
  sds = 1/(diag(out$hessian))

  if(min(eigen(out$hessian)$value)*max(eigen(out$hessian)$value)<0){
    print("bad hessian")
  }
  return(list(hessian = out$hessian
              ,parameters = name_pars
              ,sd = sds))
}




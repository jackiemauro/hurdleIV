#' Make a proper covariance matrix
#' 
#' @description In these regressions, we need to transform the covariance
#' matrix using the coefficients in our regressions. This function produces that 
#' covariance matrix.
#'
#' @param a can be a matrix or vector of what should be included in 
#' the covariance matrix
#' @param num_endog the number of endogenous regressors
#' @param gamma a vector of your second stage probit regression parameters
#' @param beta a vector of your second stage linear regression parameters
#' @param option either "mat" if a is a matrix, "vector" if a is all the
#' elements of the matrix but in vector form or "parameters" in which case
#' a should be a list of the important parameters:
#' a = list(rho,tau0,tau1,y_sd,endog_sd)
#' 
#' @return returns the covariance matrix


make.covTrans <- function(a,num_endog,gamma,beta,option = "mat",noname = F){
  #option "parameters": a = list(rho,tau0,tau1,y_sd,endog_sd)
  #option "mat": a is full matrix
  #option "vector": a is vector of all elements of matrix
  
  if(option == "mat"){
    if(is.matrix(a)){
      Sig_err = a
    }
    else{
      cat("Error: You have set option to 'mat', please supply a matrix\n
          Other options:\n
          option 'vector': Insert full matrix as a vector\n
          option 'parameters': Insert only the parameters of interest as list, unimportant off-diagonal elements will be set to 0")
      return(NA)
    }
  }
  
  else if(option=="vector"){
    tryCatch({
      Sig_err = matrix(a, ncol = 2+num_endog, byrow = F)
    }, error = function(e){
      cat("Error: You have not supplied the correct length vector \n
          Vector should be (2+number of endogenous variables)^2 long\n
          Other options:\n
          option 'mat': Insert full matrix\n
          option 'parameters': Insert only the parameters of interest as list, unimportant off-diagonal elements will be set to 0")
      return(NA)
    }, warning = function(w){
      cat("Error: You have not supplied the correct length vector \n
          Vector should be (2+number of endogenous variables)^2 long\n
          Other options:\n
          option 'mat': Insert full matrix\n
          option 'parameters': Insert only the parameters of interest as list, unimportant off-diagonal elements will be set to 0")
      return(NA)
    }
      )
    
    }
  
  else if(option == "parameters"){
    tryCatch({
      mat1 = matrix(c(1,a$rho,a$rho,a$y_sd^2),ncol = 2, byrow = T)
      tau_mat = as.matrix(cbind(a$tau0,a$tau1))
      endog_mat = diag(length(a$endog_sd))*a$endog_sd^2
      Sig_err = rbind(cbind(mat1,t(tau_mat)),cbind(tau_mat,endog_mat))
    }, error = function(e){
      print(e)
      cat("Error: You have not supplied the correct length vector \n
          Vector should include rho, y_sd, tau0, tau1, endog_sd\n
          Other options:\n
          option 'mat': Insert full matrix\n
          option 'vector': Insert full matrix as a vector\n")
      return(NA)
    }, warning = function(w){
      print(w)
      cat("Error: You have not supplied the correct length vector \n
          Vector should include rho, y_sd, tau0, tau1, endog_sd\n
          Other options:\n
          option 'mat': Insert full matrix\n
          option 'vector': Insert full matrix as a vector\n")
      return(NA)
    }
      )
    }
  
  
  A = diag(2+num_endog)
  
  if(noname == T){
    gam2 = tail(gamma,num_endog)
    bet2 = tail(beta,num_endog)
  }
  else{
    gam2 = gamma[names(gamma) %in% endog_names]
    bet2 = beta[names(beta) %in% endog_names]
  }
  
  A[1,3:(2+num_endog)] <- gam2
  A[2,3:(2+num_endog)] <- bet2
  
  Sig = A%*%Sig_err%*%t(A)
  
  if(min(eigen(Sig)$values)<=0){
    print("bad start sigma values")
    return(NA)
  }
  
  return(Sig)
  }

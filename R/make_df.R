#' Make a dataset of Normals
#'
#' @description Generates a dataset with normally distributed variables.
#' You can specify the mean and variance of each column.
#'
#' @param mean a vector of means for each column of your dataset
#' @param sd a vector of standard deviations for each column of your dataset
#' @param n the length of your variables
#' @param pref default is False but fill in if you want to specify the names
#' of your columns
#'
#' @return returns a dataset of normally distributed variables


make.df <- function(mean,sd,n, pref = F){
  if(length(mean)!=length(sd)){
    stop("Error: length of mean and sd vectors differ")
  }
  out= matrix(c(rep(NA,length(mean)*n)),ncol = length(mean))
  for(i in 1:length(mean)){
    out[,i] = rnorm(n,mean[i],sd[i])
  }

  df = as.data.frame(out)

  if(pref != F){
    names(df) <- make.names(len = length(mean),pref = pref)
  }

  return(df)
}

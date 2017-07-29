#' Helper fn
#'
#' @description Names the columns of your dataset to help out
#'
#' @param len the length of your name vector
#' @param pref the prefix
#'
#' @return returns a vector of names of the form x1, x2, .... xN


make.names <- function(len, pref){
  name = c(rep(NA,len))
  for(i in 1:len){
    name[i] = paste(pref,i,sep = "")
  }
  return(name)
}

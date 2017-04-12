#' Test if you can invert a matrix 
#'
#' @param t a matrix
#' 
#' @return True or false

cant.solve <- function(m) class(try(solve(m),silent=T))!="matrix"
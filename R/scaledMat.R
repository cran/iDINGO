scaledMat <- function(x){
 newx=x/sqrt(diag(x) %*% t(diag(x)))
 return(newx)
 }

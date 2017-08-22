trans.Fisher <-  function(x) {
  x[x>=(1-1e-07)] <- 1 - 1e-07
  x[x<=(-1+1e-07)] <- -1 + 1e-07
  return(log((1+x)/(1-x))/2)
}

extendedBIC <-
function(gamma,omegahat,S,n) {
    p = nrow(omegahat)
    es = sum(omegahat[upper.tri(omegahat)]!=0)
    return(-log(det(omegahat)) + sum(diag(omegahat %*% S)) + es * (log(n)/n) + es * gamma * (4*log(p)/n))
}

Sigmax <-
function(P=NULL,Q,Psi,x) {
    p = nrow(Psi)
    R = ncol(Q)/length(x)
    x = rep(x,R)
    II = diag(p)
    if (is.null(P)) {
        return(Q %*% x %*% t(x) %*% t(Q) + Psi) 
    }else {
        iP = solve(II-P)
        return(iP %*% (Q %*% x %*% t(x) %*% t(Q) + Psi)  %*% t(iP))}
}

single.boot <- function(i, z, n, tY.org, P, levels.z, w.upper) {
    ## i <- iteration number
    ## z <- covariate indicator (should be binary)
    ## n <- number of rows in data
    ## tY.org <- transformed standardized data
    ## P <- global component
    ## levels.z <- levels of the covariates
    ## w.upper <- upper triangular of Omega
    
w.id = sample(1:n,replace=T)
tY = tY.org[w.id,]
mdat = apply(tY,2,mean)
sdat = apply(tY,2,sd)
std.tY = t((t(tY) - mdat)/sdat)

fit.g = Greg.em(std.tY~z[w.id]) ## Fitting

smat = diag(sdat)
sigmaX1 = smat%*%Sigmax(Q=fit.g$B,P=P,Psi=fit.g$A,x=c(1,levels.z[1]))%*%smat
omegaX1 = solve(sigmaX1)
boot.RX1 = trans.Fisher(-scaledMat(omegaX1)[w.upper])

sigmaX2 = smat%*%Sigmax(Q=fit.g$B,P=P,Psi=fit.g$A,x=c(1,levels.z[2]))%*%smat
omegaX2 = solve(sigmaX2)
boot.RX2= trans.Fisher(-scaledMat(omegaX2)[w.upper])

boot.diff = trans.Fisher(-scaledMat(omegaX1)[w.upper]) - trans.Fisher(-scaledMat(omegaX2)[w.upper])

return(boot.diff)
}
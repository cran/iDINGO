scoring.boot <- function(stddat,z,Omega,A,B,boot.B=100,verbose=T) {
    ## stddat <- standardized nxp data with colnames as genename
    ## z <- covariate indicator (should be binary)
    ## Omega <- precision matrix for whole data (implies global network)
    ## A <- MLE of the baseline covariance matrix
    ## B <- MLE of the regression coefficients
    ## boot.B <- number of bootstrap samples   
    
    p = ncol(stddat)
    n = nrow(stddat)
    w.upper = which(upper.tri(Omega))
    w.mat = which(upper.tri(Omega),arr.ind=T)
    II = diag(p)
    levels.z = unique(z)
    
    stopifnot(length(levels.z)==2)
    
    # make global component from Omega
    diag.Omega = diag(Omega)
    P = -Omega/diag.Omega
    diag(P) = 0
    
    # data transformation
    tY.org = stddat %*% (II-t(P))
    
    b = 1
    boot.diff = matrix(nrow=length(w.upper),ncol=boot.B)
    while(b<=boot.B) {
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

      boot.diff[,b] = trans.Fisher(-scaledMat(omegaX1)[w.upper]) - trans.Fisher(-scaledMat(omegaX2)[w.upper])
      if (verbose & b%%10 ==0) {cat(b," bootstraps are done","\n")}
      b = b+1
    }

    # differential scoring
    genepair = data.frame(gene1 = colnames(stddat)[w.mat[,1]],gene2 =colnames(stddat)[w.mat[,2]])
    R1 = -scaledMat(solve(Sigmax(Q=B,P=P,Psi=A,x=c(1,levels.z[1]))))[w.upper]
    R2 = -scaledMat(solve(Sigmax(Q=B,P=P,Psi=A,x=c(1,levels.z[2]))))[w.upper]
    
    diff.score = (trans.Fisher(R1) - trans.Fisher(R2)) / apply(boot.diff,1,sd)
    
    # Convert differential score to corrected p-value
    p.val = getEfronp(diff.score, plotIt = FALSE)$correctp

    return(list(genepair=genepair,levels.z=levels.z,R1=R1, R2=R2, boot.diff=boot.diff,diff.score=diff.score,p.val=p.val))
}

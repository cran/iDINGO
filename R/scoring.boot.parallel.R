scoring.boot.parallel <- function(stddat,z,Omega,A,B,boot.B=100,verbose=T,cores=1) {
    ## stddat <- standardized nxp data with colnames as genename
    ## z <- covariate indicator (should be binary)
    ## Omega <- precision matrix for whole data (implies global network)
    ## A <- MLE of the baseline covariance matrix
    ## B <- MLE of the regression coefficients
    ## boot.B <- number of bootstrap samples
    ## cores <- number of cores to use 
    
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
    
    cl <- makeCluster(cores)
    
    #environment(single.boot) <- .GlobalEnv
    #environment(Greg.em) <- .GlobalEnv
    #environment(trans.Fisher) <- .GlobalEnv
    #environment(scaledMat) <- .GlobalEnv
    
    clusterExport(cl, list("single.boot", "Greg.em", "trans.Fisher", "scaledMat", "Sigmax"))
    
    # run bootstraps    
    boot.diff <-  parSapply(cl, 1:boot.B, single.boot, 
                            z=z, n=n, tY.org=tY.org, P=P, levels.z=levels.z, w.upper=w.upper)
    
    stopCluster(cl)

    # differential scoring
    genepair = data.frame(gene1 = colnames(stddat)[w.mat[,1]],gene2 =colnames(stddat)[w.mat[,2]])
    R1 = -scaledMat(solve(Sigmax(Q=B,P=P,Psi=A,x=c(1,levels.z[1]))))[w.upper]
    R2 = -scaledMat(solve(Sigmax(Q=B,P=P,Psi=A,x=c(1,levels.z[2]))))[w.upper]
    
    diff.score = (trans.Fisher(R1) - trans.Fisher(R2)) / apply(boot.diff,1,sd)
    
    # Convert differential score to corrected p-value
    p.val = getEfronp(diff.score, plotIt = FALSE)$correctp

    return(list(genepair=genepair,levels.z=levels.z,R1=R1, R2=R2, boot.diff=boot.diff,diff.score=diff.score,p.val=p.val))
}

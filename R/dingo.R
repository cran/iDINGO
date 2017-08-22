dingo <- function(dat,x,rhoarray=NULL,diff.score=T,B=100,verbose=T,cores=1) {
#  Input:
####  - dat : n by p response data
####  - x : length n covariate vector
####  - rhoarray : Candidate tuning parameters of glasso for fitting global network. If it is one value, then we use the value. If it is null, we set the candidates.
####  - diff.score : if TRUE, edge-wise differential scores are calculated from bootstrap standard error. 
####  - B : the number of bootstrap samples. It is used if diff.score is TRUE
####  - verbose : if TRUE, the number of bootstrap samples are listed
####  - cores : the number of cores to use for parallel bootstrapping
    timer <- vector(mode = "numeric", length = 3)
    ptm <- proc.time()[3]
    
    # Check that we are not specifying more cores than available. Reset if needed.
    cores.avail <- detectCores() - 1
    if (cores >= 2 & cores > cores.avail) {
      cat("Cores specified (", cores, ") is greater than recommended for your system (", cores.avail, "). Setting cores = ", cores.avail, ". \n", sep = "")
      cores <- cores.avail
    }
  
    n = nrow(dat)
    p = ncol(dat)
    II = diag(p)
    w.upper = which(upper.tri(II))
    w.mat = which(upper.tri(II), arr.ind = T)
    
    # data standardization
    mdat = apply(dat,2,mean)
    sdat = apply(dat,2,sd)
    stddat = t((t(dat) - mdat)/sdat)
    
    ###################################
    # fit global network using glasso #
    ###################################

    S = cov(stddat)
    if (is.null(rhoarray)) {
      if (n>p) {rhoarray = exp(seq(log(0.001),log(1),length=100))
      }else {rhoarray = exp(seq(log(0.1),log(3),length=100))}
    }
    BIC = rep(0,length(rhoarray))
    for (rh in 1:length(rhoarray)) {
      fit.gl1 = glasso(S,rho=rhoarray[rh])
      fit.gl2 = glasso(S,rho=rhoarray[rh],w.init=fit.gl1$w,wi.init=fit.gl1$wi)
      BIC[rh] = extendedBIC(gamma=0,omegahat=fit.gl2$wi,S=S,n=nrow(dat))
    }
   rho = rhoarray[which.min(BIC)]
   fit.gl1 = glasso(S,rho=rho)
   fit.gl2 = glasso(S,rho=rho,w.init=fit.gl1$w,wi.init=fit.gl1$wi)

   Omega = fit.gl2$wi
   diag.Omega = diag(Omega)
   P = -Omega/diag.Omega
   diag(P) = 0

   if (verbose) cat("Step 1 of DINGO is finished at",date(),"\n")
   timer[1] <- proc.time()[3] - ptm

   ######################################
   # fit local group-specific component #
   ######################################
   tY.org = stddat %*% (II-t(P))
   mdat = apply(tY.org,2,mean)
   sdat = apply(tY.org,2,sd)
   std.tY = t((t(tY.org) - mdat)/sdat)
   fit.g = Greg.em(std.tY~x)
   
   if (verbose) cat("Step 2 of DINGO is finished at",date(),"\n")
   timer[2] <- proc.time()[3] - timer[1] - ptm
   
   if (diff.score) {
     if (verbose) cat("Bootstrap scoring is started at", date(),"\n")
     if (cores >= 2) {
       boot.fit = scoring.boot.parallel(stddat=stddat,z=x,Omega=Omega,A=fit.g$A,B=fit.g$B,boot.B=B,verbose=verbose,cores=cores)
     } else {
       boot.fit = scoring.boot(stddat=stddat,z=x,Omega=Omega,A=fit.g$A,B=fit.g$B,boot.B=B,verbose=verbose)
     }
     if (verbose) cat("Bootstrap scoring is done at", date(),"\n")
     timer[3] <- proc.time()[3] - timer[2] - ptm
     
     return(list(genepair = boot.fit$genepair,levels.x = boot.fit$levels.z,R1=boot.fit$R1,R2=boot.fit$R2,boot.diff=boot.fit$boot.diff,diff.score=boot.fit$diff.score,p.val=boot.fit$p.val,rho=rho,P=P,Q=fit.g$B,Psi=fit.g$A,step.times=timer))
   }else{
     genepair = data.frame(gene1 = colnames(dat)[w.mat[, 1]], 
        gene2 = colnames(stddat)[w.mat[, 2]])
     levels.x = unique(x)
     R1 = -scaledMat(solve(Sigmax(Q = fit.g$B, P = P, Psi = fit.g$A, x = c(1,
        levels.x[1]))))[w.upper]
     R2 = -scaledMat(solve(Sigmax(Q = fit.g$B, P = P, Psi = fit.g$A, x = c(1, 
        levels.x[2]))))[w.upper]
     return(list(genepair = genepair,levels.x = levels.x,R1=R1,R2=R2,boot.diff=NULL,diff.score=NULL,p.val=NULL,rho=rho,P=P,Q=fit.g$B,Psi=fit.g$A,step.times=timer))
   }
}


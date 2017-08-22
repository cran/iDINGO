Greg.em <-
function(formula,data=NULL,R=1,tol=1e-10,itmax=1000,verbose=F){
    model=model.frame(formula,data)
    Y=model.response(model)
    X=model.matrix(formula,model)
    
    ## normal log likelihood ##
    
    ldmvnorm=function(y,mu=rep(0,length(y)),Sig=diag(1,length(y))){
        -.5*( length(y)*log(2*pi) + log(det(Sig)) + t(y-mu)%*%solve(Sig)%*%(y-mu)  )
    }
    
    p=dim(Y)[2] ; q=dim(X)[2] ; n=dim(Y)[1]
    
    # Initial values
    dvar = apply(Y,2,var)
	A = diag(dvar)
	iA=diag(1/dvar)
	B=matrix(rep(t(coef(lm(Y~X[,2]))),R),ncol=q*R)
    
    # Fitting
	iter=0
	LL=NULL
	rll=10
    
    while( rll > tol  & iter<itmax)
    {
        B0=B ; iter=iter+1
        
        ### find expectation, var of z
        Vz=array(dim=c(R,R,n)) ; Mz=matrix(nrow=n,ncol=R)
        for(i in 1:n)
        {
            Bx=apply(array(B,dim=c(p,q,R)),3,"%*%",X[i,])
            Vz[,,i]=solve(  t(Bx)%*%iA%*%Bx + diag(R) )
            Mz[i,]=Vz[,,i]%*%t(Bx)%*%iA%*%Y[i,]
        }
        ###
        
        ### obtain MLEs
        Y1=Y ; X1=NULL ; for(r in 1:R) { X1=cbind(X1,diag(Mz[,r])%*%X  )}
        Y0=matrix(0,nrow=n*R,ncol=p) ; X0=NULL
        for(i in 1:n)
        {
            xi=matrix(outer(X[i,],diag(R)),nrow=R*q,ncol=R)
            ZZ=xi%*%Vz[,,i]%*%t(xi) ; ZZ=.5*(ZZ+t(ZZ))
            Z=eigen(ZZ);Z=Z$vec[,1:R]%*%diag(sqrt(Z$val[1:R]),nrow=R)
            X0=rbind(X0,t(Z))
        }
        YA=rbind(Y0,Y1) ; XA=rbind(X0,X1)
        
        B=t(YA)%*%XA%*%solve(t(XA)%*%XA)
        E=YA-XA%*%t(B)
        A= (t(E)%*%E)/n
        dA=diag(A)
        iA=diag(1/dA)
        A = diag(dA)
        ###
        
        
        ###
        if(iter%%5==0)
        {
            ll=0
            for(i in 1:dim(Y)[1])
            {
                xi=matrix(outer(X[i,],diag(R)),nrow=R*q,ncol=R)
                ll=ll+ldmvnorm(Y[i,],Sig=A+B%*%xi%*%t(xi)%*%t(B))
            }
            LL=c(LL,ll)
            if(iter>5){rll=abs(LL[length(LL)]-LL[length(LL)-1])/abs(LL[length(LL)])}
            if (verbose )cat(iter,log(rll,base=10),ll," ",round(diag(A),2)," ",round(c(B),2),"\n")
        }
        ###
    }
    list(A=A,B=B)
}

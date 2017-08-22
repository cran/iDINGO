\name{Sigmax}
\alias{Sigmax}
\title{
 group specific covariance matrices
}
\description{
From parameters of DINGO model, group specific covariance matrices are obtained
}
\usage{
Sigmax(P = NULL, Q, Psi, x)
}
\arguments{
  \item{P}{
   a p x p matrix specifying global component
}
  \item{Q}{
   the coefficient parameter matrix of covariance regression model using Greg.em function
}
  \item{Psi}{
   the diagonal error variance matrix of covariance regression model using Greg.em function
}
  \item{x}{
   a vector specifying group. This must be corresponding to the design matrix of Greg.em function
}
}
\value{
group specific precision matrix}
\author{
Min Jin Ha <mjha@mdanderson.org>
}
\examples{
library(glasso)
data(gbm)
x = gbm[,1]
Y = as.matrix(gbm[,-1])
p = ncol(Y)
# Estimating inverse covariance matrix using GLasso #
S = cov(Y)
w.upper = which(upper.tri(S))

rhoarray = exp(seq(log(0.001),log(1),length=100))
BIC = rep(0,length(rhoarray))
for (rh in 1:length(rhoarray)) {
    fit.gl1 = glasso(S,rho=rhoarray[rh])
    BIC[rh] = extendedBIC(gamma=0,omegahat=fit.gl1$wi,S=S,n=nrow(Y))
}
rho = rhoarray[which.min(BIC)]
fit.gl2 = glasso(S,rho=rho)
Omega = fit.gl2$wi

# Fitting (Covariance Regression on transformed data)
diag.Omega = diag(Omega)
P = -Omega/diag.Omega
diag(P) = 0

tY = Y %*% (diag(p)-t(P))
mdat = apply(tY,2,mean)
sdat = apply(tY,2,sd)
std.tY = t((t(tY) - mdat)/sdat)
smat = diag(sdat)

## rank 1 covariance regression
fit.g = Greg.em(std.tY~x,R=1) 
## obtain covariance matrix of Y when x=1
sigmaX1 = Sigmax(Q=fit.g$B,P=P,Psi=fit.g$A,x=c(1,1)) 
}

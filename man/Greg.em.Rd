\name{Greg.em}
\alias{Greg.em}
\title{Fitting precision regression models
}
\description{
This function fits the covariance regression model by Hoff and Niu (2012) using EM algorithm with the restriction of diagonal matrix for the noise variance
}
\usage{
Greg.em(formula, data = NULL, R = 1, tol = 1e-10, itmax = 1000, verbose = F)
}
\arguments{
  \item{formula}{
an object of class "formula" used in model.frame function}
  \item{data}{
a data frame used in model.frame function
}
  \item{R}{
  rank of the model 
}
  \item{tol}{
 a stopping criterion
}
  \item{itmax}{
  maximum number of iteration
}
  \item{verbose}{
   If true, estimation results for each iteration are printed
}
}
\value{
  \item{A}{MLE of the baseline covariance matrix}
  \item{B}{MLE of the regression coefficients}
}
\references{
Hoff, P. D. and Niu, X. (2012) A covariance regression model. Statistica Sinica, 22, 729-753.
}
\author{
Min Jin Ha <mjha@mdanderson.org>}
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
}


\name{extendedBIC}
\alias{extendedBIC}
\title{Extended bayesian information criteria for gaussian graphical models
}
\description{
Extended bayesian information criteria for gaussian graphical models
}
\usage{
extendedBIC(gamma,omegahat,S,n)
}
\arguments{
  \item{gamma}{
    a tuning parameter taking a scalar in [0,1] and leading to stronger penalization of large graphs
}
  \item{omegahat}{
    a p x p matrix indicating an estimates of precision (inverse covariance) matrix
}
  \item{S}{
    a p x p matrix indicating sample covariance matrix
}
 \item{n}{
 a scalar indicating sample size   
}
}
\value{
    Extended BIC penalized by the size of graphs
}
\references{
   Foygel, R. and Drton, M. (2010). Extended bayesian information criteria for gaussian graphical models. arXiv preprint arXiv:1011.6640 .
}
\author{
Min Jin Ha <mjha@mdanderson.org>}
\examples{
library(glasso)
data(gbm)
x = gbm[,1]
Y = gbm[,-1]

# Estimating inverse covariance matrix using GLasso #
S = cov(Y)

rhoarray = exp(seq(log(0.001),log(1),length=100))
BIC = rep(0,length(rhoarray))
for (rh in 1:length(rhoarray)) {
    fit.gl1 = glasso(S,rho=rhoarray[rh])
    BIC[rh] = extendedBIC(gamma=0,omegahat=fit.gl1$wi,S=S,n=nrow(Y))
}
rho = rhoarray[which.min(BIC)]
fit.gl2 = glasso(S,rho=rho)
Omega = fit.gl2$wi


}

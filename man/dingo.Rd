\name{dingo}
\alias{dingo}
\title{
Fit DINGO model
}
\description{
This function fit DINGO model and calculates edge-wise differential scores for all pairwise edges among p variables. 
}
\usage{
dingo(dat,x,rhoarray=NULL,diff.score=T,B=100,verbose=T,cores=1)
}
\arguments{
  \item{dat}{
       nxp data with colnames as genename
  }
  \item{x}{a length n vector representing a binary covariate
  }
  \item{rhoarray}{a vector representing candidate tuning parameters of glasso for fitting global network model. If it is one value, then we use the value as the tuning parameter. It is set by NULL as default and we select 100 candidate values.
  }
  \item{diff.score}{a logical value. If TRUE, edge-wise differential scores are calculated from bootstrap standard error. Otherwise, we fit Steps 1 and 2 of DINGO model to get group specific GGMs (partial correlations)  }
  \item{B}{the number of bootstrap samples to calculate differential scores.
  }
  \item{verbose}{if TRUE, lists the procedure
  }
  \item{cores}{the number of cores to run in parallel for bootstrapping, set to 1 as a default. If more cores are specified than the recommended maximum (the number of cores detected minus 1), this value will be replaced by the recommended value. 
  }
}
\value{
    \item{genepair}{a p(p-1)/2 x 2 matrix indicating all pairs of genes
    }
    \item{levels.x}{a length 2 vector indicating levels of the binary covariate x, the first element is for group 1 and the second element is for group 2
    }
    \item{R1}{a length p(p-1)/2 vector indicating partial correlations for group 1 and the order is corresponding to the order of genepair
    }
    \item{R2}{a length p(p-1)/2 vector indicating partial correlations for group 2 and the order is corresponding to the order of genepair
    }
    \item{boot.diff}{a p(p-1)/2 x boot.B matrix indicating bootstrapped difference, Fisher's Z transformed R1 - R2. The rows are corresponding to the order of gene pair and the columns are corresponding to the bootstrap samples
    }
    \item{diff.score}{a p(p-1)/2 vector of differential score corresponding to genepair
    }
    \item{p.val}{a p(p-1)/2 vector of corrected p-values corresponding to genepair
    }
    \item{rho}{selected tuning parameter of glasso fit
    }
    \item{P}{p by p matrix of Global component of the DINGO model
    }
    \item{Q}{p by 2 matrix of the coefficient parameter of the local group specific component L(x) of the DINGO model.}
    \item{Psi}{ p by p diagonal matrix of the noise covariance parameter of the local group specific component L(x) of the DINGO model.
    }
    \item{step.times}{ a length 3 vector containing the elapsed time for Step 1, Step 2, and Bootstrap Scoring, respectively.
    }
}
\author{
    Min Jin HA mjha@mdanderson.org, Caleb CLASS caclass@mdanderson.org	
}


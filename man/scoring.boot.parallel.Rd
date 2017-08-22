\name{scoring.boot.parallel}
\alias{scoring.boot.parallel}
\title{
Calculating differential score with parallel bootstrap scoring
}
\description{
This function calculates standard errors for edge-wise partial correlation differences obtained from DINGO model. Bootstrapping is done in parallel using parSapply from the "parallel" library.
}
\usage{
scoring.boot.parallel(stddat,z,Omega,A,B,boot.B=100,verbose=T,cores=1)
}
\arguments{
  \item{stddat}{
       standardized nxp data with colnames as genename
  }
  \item{z}{a length n vector representing a binary covariate
  }
  \item{Omega}{a p x p precision matrix for std dat which implies the global network
  }
  \item{A}{p x p matrix of the MLE for the baseline covariance matrix which is obtained from A value of the Greg.em function.
  }
  \item{B}{p x 2 matrix of the MLE for the regression coefficient which is obtained from B value of the Greg.em function
  }
  \item{boot.B}{a scalar indicating the number of bootstraps
  }
  \item{verbose}{if TRUE, lists the bootstrap replications
  }
  \item{cores}{the number of cores to run in parallel for bootstrapping, set to 1 as a default.
  }
}
\value{
    \item{genepair}{a p(p-1)/2 x 2 matrix indicating all pairs of genes
    }
    \item{levels.z}{a length 2 vector indicating levels of the binary covariate z, the first element is for group 1 and the second element is for group 2
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
}
\author{
    Min Jin HA mjha@mdanderson.org, Caleb CLASS caclass@mdanderson.org
}


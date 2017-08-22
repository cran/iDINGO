\name{idingo}
\alias{idingo}
\title{
Fit iDINGO model
}
\description{
This function fits the iDINGO model and calculates edge-wise differential scores for all pairwise edges among p variables between multiple platforms.
}
\usage{
idingo(dat,dat2=NULL,dat3=NULL,x,plats=NULL,rhoarray=NULL,
        diff.score=T,B=100,verbose=T,cores=1)
}
\arguments{
  \item{dat}{
       nxp dataframe/matrix with colnames as genename
  }
  \item{dat2}{
       Second nxp dataframe/matrix with colnames as genename (optional)
  }
  \item{dat3}{
       Third nxp dataframe/matrix with colnames as genename (optional)
  }
  \item{x}{a length n vector representing a binary covariate
  }
  \item{plats}{a length 1-3 vector (corresponding to the number of data sets submitted, with names for the platforms/levels of the data, such as "microRNA" or "RNAseq". This is optional, and default names "platN" will be used if names are not provided.
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
    \item{diff.score}{a p(p-1)/2 vector of differential score corresponding to genepair
    }
    \item{p.val}{a p(p-1)/2 vector of corrected p-values corresponding to genepair
    }
}
\author{
    Caleb CLASS caclass@mdanderson.org, Min Jin HA mjha@mdanderson.org
}


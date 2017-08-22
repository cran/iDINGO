\name{single.boot}
\alias{single.boot}
\title{
Calculating differential score for a single bootstrap
}
\description{
This function calculates the edge-wise partial correlation difference for a single bootstrap.
}
\usage{
single.boot(i, z, n, tY.org, P, levels.z, w.upper)
}
\arguments{
  \item{i}{iteration number. This is not used within this function, but necessary for parSapply within scoring.boot.parallel function.
  }
  \item{z}{a length n vector representing a binary covariate
  }
  \item{n}{the number of rows in data
  }
  \item{tY.org}{the transformed standardized data
  }
  \item{P}{the global correlation component
  }
  \item{levels.z}{the levels of the covariates
  }
  \item{w.upper}{the upper triangular of Omega
  }
}
\value{
    \item{boot.diff}{the difference for this bootstrap
    }
}
\author{
    Min Jin HA mjha@mdanderson.org, Caleb CLASS caclass@mdanderson.org
}



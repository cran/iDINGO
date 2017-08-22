\name{plotNetwork}
\alias{plotNetwork}
\title{
Plot differential network
}
\description{
This function plots the differential network from a completed DINGO or iDINGO model.
}
\usage{
plotNetwork(fit, threshold=0.05, thresh.type="p.val", layout="circular", 
            legend.pos="left")
}
\arguments{
  \item{fit}{
       output from running dingo() or idingo()
  }
  \item{threshold}{
       a numeric value containing the threshold for which edges will be included in the differential network plot. If 'thresh.type' is 'p.val', all edges with p-values below this threshold will be included in the plot. If 'thresh.type' is 'diff.score', all edges with absolute differential scores above this threshold will be included in the plot.
  }
  \item{thresh.type}{either 'p.val' or 'diff.score', defining which variable is used as threshold for edge inclusion.
  }
  \item{layout}{either 'circular' or one of igraph's supported layouts. If 'circular', dingo networks will be plotted in a circle, and idingo networks will be plotted as a cylinder (with each platform/level as a separate circle).
  }
  \item{legend.pos}{Legend position for multi-platform networks, in c("left", "right"). Legend is not included for single-platform networks.}
}
\value{
    \item{visNet}{a network plot, using igraph and visNetwork
    }
}
\note{
To calculate differential scores and p-values for use in network plot thresholding, diff.score must be set to TRUE in dingo() or idingo().
}
\author{
    Caleb CLASS caclass@mdanderson.org, Min Jin HA mjha@mdanderson.org
}


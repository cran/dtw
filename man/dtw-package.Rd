\name{dtw-package}
\alias{dtw-package}

\docType{package}
\title{
Dynamic Time Warp algorithms in R
}
\description{
  Dynamic Time Warp: find the optimal alignment between two time series.
}

\details{
\tabular{ll}{
Package: \tab dtw\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2007-12-06\cr
License: \tab GPL-2\cr
}

Comprehensive implementation of Dynamic Time Warping algorithms in R.
Most variants of the algorithm are supported: symmetric, asymmetric and
custom step patterns, with weighting (see \code{\link{stepPattern}}).
Supports windowing: none, "Itakura" parallelogram, Sakoe-Chiba band,
custom (see \code{\link{dtwWindowingFunctions}}).  Handles query and
template of arbitrary lengths. Multivariate matching and arbitrary
definition for a distance function are supported via user-supplied local
distance matrix.

Currently pure-R.

Please see documentation for function \code{\link{dtw}}.

You are encouraged to download the latest version at
\url{http://r-forge.r-project.org/projects/dtw/}. The command
\code{install.packages("dtw",repos="http://r-forge.r-project.org")}
should get you the latest version.

}
\author{
  Toni Giorgino, Copyright (c) 2007 

  Maintainer:  toni.giorgino@unipv.it
}

\seealso{  \code{\link[pkg:analogue]{distance}},
  \code{\link{dtwWindowingFunctions}},  \code{\link{stepPattern}} , \code{\link{outer}} }

\references{
   TODO
}

\examples{
 library(dtw);
 ## demo(dtw);
}


\keyword{ package }
\keyword{ ts }

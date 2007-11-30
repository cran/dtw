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
Version: \tab 0.4\cr
Date: \tab 2007-12-10\cr
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

Please see documentation for function \code{\link{dtw}}, which is the
main entry point to the package.

A fast, compiled version of the algorithm is normally used.
Should it not be available, a slower pure-R equivalent is 
automatically used as a fall-back.



The package home page is at \url{http://dtw.r-forge.r-project.org}. The
command
\code{install.packages("dtw",repos="http://r-forge.r-project.org")}
should get you the development (possibly unstable) version.

}
\author{
  Toni Giorgino, Copyright (c) 2007 

  Maintainer:  toni.giorgino@unipv.it
}

\seealso{ \code{\link{dtw}} for the main entry point to the package;
  \code{\link{dtwWindowingFunctions}} for global constraints;
  \code{\link{stepPattern}} for local constraints;
  \code{\link[pkg:analogue]{distance}},  \code{\link{outer}} for
  building a local cost matrix with multivariate
  timeseries and custom distance functions. }

\references{
   TODO
}

\examples{
 library(dtw);
 ## demo(dtw);
}


\keyword{ package }
\keyword{ ts }

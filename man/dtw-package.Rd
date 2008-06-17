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
Version: \tab 1.0\cr
Date: \tab 2007-12-10\cr
License: \tab GPL-2\cr
}

Comprehensive implementation of Dynamic Time Warping (DTW) algorithms in
R.

DTW finds the optimal (least cumulative distance) mapping between a
given query into a given template time series.

Most variants of the algorithm are supported: symmetric, asymmetric and
custom step patterns, with weighting (see \code{\link{stepPattern}}).
Supports windowing: none, "Itakura" parallelogram, Sakoe-Chiba band,
custom (see \code{\link{dtwWindowingFunctions}}).  Handles query and
template of arbitrary lengths. Multivariate matching and arbitrary
definition for a distance function are supported via user-supplied local
distance matrix. The Minimum Variance Matching algorithm is also
supported, as a special case of DTW.

Package provides minimum cumulative distance, warping function, plots,
etc. A fast, compiled version of the algorithm is normally used.  Should
it not be available, a slower pure-R equivalent is automatically used as
a fall-back.


Please see documentation for function \code{\link{dtw}}, which is the
main entry point to the package.

If you use this software, please cite it according to
\code{citation("dtw")}.  The package home page is at
\url{http://dtw.r-forge.r-project.org}.

To get the latest stable version from CRAN, use
\code{install.packages("dtw")}.  To get the development version
(possibly unstable), use
\code{install.packages("dtw",repos="http://r-forge.r-project.org")}.

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

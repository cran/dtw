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
Date: \tab 2007-11-30\cr
License: \tab GPL-2\cr
}

Implementation of Dynamic Time Warping algorithms in R.  Most variants
of the algorithm are supported: currently symmetric, asymmetric and
custom step patterns, with weighting. Also supports windowing (Itakura
parallelogram, Sakoe-Chiba band, custom, no windowing), and handles
cases when query and template length differ.

Multivariate matching supported via user-supplied distance matrix.

Currently pure-R.

Please see documentation for function \code{\link{dtw}}.

}
\author{
  Toni Giorgino

  Maintainer:  toni.giorgino@unipv.it
}
\references{
   TODO
}
\keyword{ package }
\keyword{ ts }
\examples{
 ## TODO
}

\name{dtwDist}
\alias{dtwDist}
\title{Compute a dissimilarity matrix}

\description{Compute the dissimilarity matrix  between
  a set of single-variate timeseries.  }

\usage{
dtwDist(m,...)
# dist(m,y,method="DTW",...)
}
\arguments{
  \item{m}{numeric matrix, containing timeseries as rows}
  \item{y}{numeric matrix, containing timeseries as rows (for cross-distance)}
  \item{...}{arguments passed to the \code{\link{dtw}} call}
}
\value{
  A square matrix whose element \code{[i,j]} holds the Dynamic Time Warp
 distance between row \code{i} (query) and \code{j} (template) of
 \code{m}, i.e.  \code{dtw(m[i,],m[j,])$distance}.
}
\details{

  The \code{dtwDist} command is obsolete and has been superseded by the
  \code{\link[pkg:proxy]{dist}} function of package \pkg{proxy}; the DTW
  distance is registered as \code{method="DTW"} (see examples below).

  For asymmetric variants, make a \code{crossdist} object with the
  two-arguments version of \code{dist}.

  \code{dtwDist} computes a dissimilarity matrix, akin to
  \code{\link{dist}}, based on the Dynamic Time Warping definition of a
  distance between single-variate timeseries.

  The function returns a
  square matrix, whereas the \code{dist} object is
  lower-triangular. This makes sense because in general the DTW
  "distance" is not symmetric (see e.g.  asymmetric step patterns).
  If a proper \code{\link{dist}} object is desired, a suitable
  conversion strategy has to be chosen (see examples).
 
}


\seealso{Other "distance" functions are: \code{\link{dist}},
  \code{\link[pkg:vegan]{vegdist}} in package \code{vegan},
  \code{\link[pkg:analogue]{distance}} in package \code{analogue}, etc.
}

\examples{

## Symmetric step pattern => symmetric dissimilarity matrix;
## no problem coercing it to a dist object:

m <- matrix(0,ncol=3,nrow=4)
m <- row(m)
dist(m,method="DTW");

# Old-fashioned call style would be:
#   dtwDist(m)
#   as.dist(dtwDist(m))



## Asymmetric step pattern: we can either disregard part of the pairs
## (as.dist), or average with the transpose

mm <- matrix(runif(12),ncol=3)
dm <- dist(mm,mm,method="DTW",step="asymmetric"); # a crossdist object

# Old-fashioned call style would be:
#   dm <- dtwDist(mm,step=asymmetric)
#   as.dist(dm)


## Symmetrize by averaging:
(dm+t(dm))/2


## check definition
stopifnot(dm[2,1]==dtw(mm[2,],mm[1,],step=asymmetric)$distance)


}

\author{Toni Giorgino}
\keyword{ts}
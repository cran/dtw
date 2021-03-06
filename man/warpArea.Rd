% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/warpArea.R
\name{warpArea}
\alias{warpArea}
\title{Compute Warping Path Area}
\usage{
warpArea(d)
}
\arguments{
\item{d}{an object of class \code{dtw}}
}
\value{
The area, not normalized by path length or else.
}
\description{
Compute the area between the warping function and the diagonal (no-warping)
path, in unit steps.
}
\details{
Above- and below- diagonal unit areas all count \emph{plus} one (they do not
cancel with each other).  The "diagonal" goes from one corner to the other
of the possibly rectangular cost matrix, therefore having a slope of
\code{M/N}, not 1, as in \code{\link[=slantedBandWindow]{slantedBandWindow()}}.

The computation is approximate: points having multiple correspondences are
averaged, and points without a match are interpolated. Therefore, the area
can be fractionary.
}
\note{
There could be alternative definitions to the area, including
considering the envelope of the path.
}
\examples{

  ds<-dtw(1:4,1:8);

  plot(ds);lines(seq(1,8,len=4),col="red");

  warpArea(ds)

  ## Result: 6
  ##  index 2 is 2 while diag is 3.3  (+1.3)
  ##        3    3               5.7  (+2.7)
  ##        4   4:8 (avg to 6)    8   (+2  )
  ##                                 --------
  ##                                     6

}
\author{
Toni Giorgino
}
\keyword{ts}

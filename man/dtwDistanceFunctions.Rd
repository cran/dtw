\name{dtwDistanceFunctions}
\alias{dtwDistanceFunctions}
\alias{euclideanSquared}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distance functions for DTW}
\description{
  Compute Euclidean squared distance between two given vectors
  (default function for \code{dtw}).
}
\usage{
euclideanSquared(a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{ The first element }
  \item{b}{ The second element }
}
\details{
  This is a trivial distance function supplied as an example of what can be
  passed to the \code{distance.function} argument in \code{\link{dtw}}.
  
  Altough it makes  sense for \code{a} and \code{b} to be
  \emph{vectors} rather than numbers, direct alignment of two
  multivariate time series is currently not supported natively by
  \code{dtw} (this may change soon).  The user should instead build a local
  distance matrix, e.g. via a custom function or
  \code{\link[pkg:analogue]{distance}} in package \code{analogue}, and
  feed the result into \code{dtw} as a matrix first argument.

}

\value{
  Local distance (float).
}
\references{Euclid's Elements, T.L. Heath (translator), Dover (1956).}
\author{Toni Giorgino }


\seealso{  See Also \code{\link[pkg:analogue]{distance}}  }
\examples{

# A constant vector
z<-numeric(10);     

# dtw between two constant vectors with a squared euclidean distance
# this should be 2^2 times 10 elements = 40
d2<-dtw(z,z+2);
stopifnot(d2$distance==40);

# dtw between two constant vectors with a root squared euclidean distance
# this should be 2 times 10 elements = 20
d1<-dtw(z,z+2,
 	distance.function=function(a,b){return(abs(a-b));})
stopifnot(d1$distance==20);


}


\keyword{ internal }

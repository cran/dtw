\name{euclideanSquared}
\alias{euclideanSquared}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Euclidean Distance between two vectors }
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
  
  Altough it makes a lot of sense for \code{a} and \code{b} to be
  \emph{vectors} rather than numbers, direct alignment of two
  multivariate time series is currently not supported natively by
  \code{dtw} (this may change soon).  The user can instead build a local
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
## The function is currently defined as
function(a,b) {
  z<-a-b;
  z<-drop(z \%*\% z);                     #inner dot product
  return (z);
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }

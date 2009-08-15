\name{mvm}
\alias{mvm}
\alias{mvmStepPattern}

\title{Minimum Variance Matching algorithm}

\description{Step patterns to compute the Minimum Variance Matching
  (MVM) correspondence between  time series}

\usage{
  mvmStepPattern(elasticity=20);
}
\arguments{
  \item{elasticity}{integer: maximum consecutive reference elements skippable}
}
\value{
  A step pattern object.
}

\details{

  The Minimum Variance Matching algorithm [1] finds the non-contiguous
  parts of reference which best match the query, allowing for arbitrarily
  long "stretches" of reference to be excluded from the match. All
  elements of the query have to be matched. First and last elements of
  the query are anchored at the boundaries of the reference.

  The \code{mvmStepPattern} function creates a \code{stepPattern} object
  which implements this behavior, to be used with the usual
  \code{\link{dtw}} call (see example). MVM is computed as a special
  case of DTW, with a very large, asymmetric-like step pattern.

  The \code{elasticity} argument limits the maximum run length of
  reference which can be skipped at once. If no limit is desired, set
  \code{elasticity} to an integer at least as large as the reference
  (computation time grows linearly).
 
}


\seealso{Other objects in \code{\link{stepPattern}}.}


\references{
[1] Latecki, L. J.; Megalooikonomou, V.; Wang, Q. & Yu, D. \emph{An elastic
partial shape matching technique} Pattern Recognition, 2007, 40,
3069-3080 
\cr \cr
[2]  Toni Giorgino. \emph{Computing and Visualizing Dynamic Time Warping
  Alignments in R: The dtw Package.}  Journal of Statistical
  Software, 31(7), 1-24. \url{http://www.jstatsoft.org/v31/i07/}
  }

\examples{

## The hand-checkable example given in ref. [1] above
diffmx <- matrix( byrow=TRUE, nrow=5, c(
  0,  1,  8,  2,  2,  4,  8,
  1,  0,  7,  1,  1,  3,  7,
 -7, -6,  1, -5, -5, -3,  1,
 -5, -4,  3, -3, -3, -1,  3,
 -7, -6,  1, -5, -5, -3,  1 ) ) ;

## Cost matrix
costmx <- diffmx^2;

## Compute the alignment
al <- dtw(costmx,step.pattern=mvmStepPattern(10))

## Elements 4,5 are skipped
print(al$index2)

plot(al,main="Minimum Variance Matching alignment")



}

\author{Toni Giorgino}
\keyword{ts}
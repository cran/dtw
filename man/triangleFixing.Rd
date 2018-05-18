\name{triangleFixing}
\alias{triangleFixing}
\alias{tri.ineq.show}

\title{Triangle Fixing algorithm for Metric Nearness}
\description{
  Return the matrix satisfying the triangle inequality
  most similar to a given input matrix.
}
\usage{
triangleFixing(D, tolerance=1e-15, max.iter=1e6)
tri.ineq.show(D)
}
\arguments{
  \item{D}{the matrix to be fixed}
  \item{tolerance}{threshold for convergence}
  \item{max.iter}{iteration limit}
}

\value{
  The fixed distance matrix (square symmetric).
}
\details{ Given a symmetric matrix, returns 
  the symmetric matrix "most similar" to it that satisfies
  the triangle inequality. (Similarity is measured in the L2 sense.)
  This is an implementation of the L2 triangle fixing algorithm
  for metric nearness (\code{Metric_Nearness_L2}) given in [1] and [2].

  WARNING: the algorithm does not seem to work in all cases.

  The input should be a symmetric matrix (or dist object, which is
  coerced). Convergence of the  iterative algorithm is assumed 
  when total matrix element change is less than \code{tolerance} or the
  \code{max.iter} iteration limit is reached.

  The \code{tri.ineq.show} function checks whether a dissimilarity
  matrix respects the triangle inequality. Returns NULL if it does, or
  a Nx3 matrix of indices where it doesn't (unsorted).
}
\references{
  1. Brickell, J., Dhillon, I., Sra, S., Tropp, J. (2008). The Metric Nearness Problem. SIAM. J. Matrix Anal. & Appl. 30, 375-396.

  2. S. Sra, "METRIC: The Metric Nearness Problem" code at \url{http://suvrit.de/work/soft/metricn.html}.
}
\seealso{
\code{\link[fossil]{tri.ineq}} in package \code{fossil}, which checks
whether a distance matrix respects the triangle inequality.
}

\examples{
%% dyn.load("src/dtw.so"); source("R/triangleFixing.R"); 
a<-matrix(0.2, 4,4)
a[4,2]<-0.8
a<-as.matrix(as.dist(a))
af<-triangleFixing(a)

tri.ineq.show(a)
stopifnot(tri.ineq.show(af)==NULL)

## Example in http://suvrit.de/work/soft/metricn.html
a<-matrix(c(0, 3, 7,
            3, 0, 1,
            7, 1, 0),3)
af<-triangleFixing(a)

tri.ineq.show(a)
stopifnot(tri.ineq.show(af)==NULL)

}

\concept{Dissimilarity matrix}
\concept{Triangle inequality}

\author{Toni Giorgino}
\keyword{cluster}

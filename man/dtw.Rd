\name{dtw}
\alias{dtw}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Dynamic Time Warp}
\description{
  Compute Dynamic Time Warp
  and find optimal alignment between two time series.
}
\usage{
dtw(x, y=NULL, distance.function=euclideanSquared, step.pattern="s",
  window.type="none", window.size=10, keep.internals=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Query vector OR local cost matrix }
  \item{y}{ Template vector, or unused if cost matrix given }
  %  \item{partial}{ ~~Describe \code{partial} here~~ }
  \item{distance.function}{ Pointwise distance function. See e.g. \code{\link{euclideanSquared}} }
  \item{step.pattern}{ The step pattern. Char: "s" (symmetric), "a"
    (asymmetric), or an $s\times3$ matrix containing the allowed
    steps with their cost -- to be documented better}
  \item{window.type}{ Windowing function. Char: "none" (none, default),
    "itakura" (Itakura Parallelogram), "sakoechiba" (Sakoe-Chiba band),
    or a  function}
  \item{window.size}{Window size (used by some
    \cite{window.type}s)}
  \item{keep.internals}{Don't discard the cumulative cost matrix and other
    internal structures}
}

\details{

  The function performs Dynamic Time Warp (DTW) and computes the optimal
  alignment between two time series \code{x} and \code{y}, given as
  numeric vectors.  The ``optimal'' alignment minimizes the sum of
  distances between aligned elements.
  
  The (local) distance between elements of \code{x} (query) and \code{y}
  (template) is computed through the \code{distance.function}
  (see). Lengths of \code{x} and \code{y} may  differ.

  A local-distance matrix can be passed, as well; element $i,j$ of the
  local-distance matrix is the distance between element x[i] and
  y[j]. The distance matrix has therefore $n=length(x)$ rows and
  $m=length(y)$ columns (see note below).

  Given the two time series, a local distance matrix can also be
  computed by the user, for example e.g. via the function \code{outer}
  or the complete \code{\link[pkg:analogue]{distance}} in package
  \code{analogue}.

  Common variants of DTW are supported via the \code{step.pattern}
  argument, namely: 

  \preformatted{
    DTW Type     step.pattern     step.pattern (as matrix)
    --------------------------------------------------
      Symmetric    "s"              [1,] c(1,0,1)   # advance in x
                                    [2,] c(0,1,1)   # advance in y
                                    [3,] c(1,1,1)   # diagonal

      Asymmetric   "a"              [1,] c(1,0,1)   # duplicate
                                    [2,] c(1,1,1)   # copy
				    [3,] c(1,2,1)   # skip
  }

  A custom step pattern is a matrix with three columns. Each row has an
  allowed step size in the $x$ (element [,1]) and $y$ (element [,2])
  series. Element [,3] is the weight which multiplies the local
  distance, if that step is taken.  Steps must be \emph{non-negative}
  integers (otherwise the monotonicity constraint is violated); weight
  can be float.
				  

  Windowing is supported by supplying a name into the \code{window.type}
  argument (abbreviations allowed)
  between:
  \preformatted{
    "none"         No windowing (default)
    "sakoechiba"   A "diagonal" band (see below if rectangular)
    "itakura"      Itakura parallelogram (unimplemented yet)
  }

  Alternatively, \code{window.type} can be a boolean function of two
  integer arguments to return \code{TRUE} if the given coordinates fall
  within the window (see e.g.\    \code{\link{sakoeChibaWindow}}).
  
  If the cost matrix is non-square, \emph{built-in} window types are
  automatically adjusted. The Sakoe-Chiba band, for example, is centered
  about the (jagged) line which joins element $(1,1)$ to element $(n,m)$,
  and \code{window.size} columns wide.


  

  TODO references
  

  
}
\value{
  An object (currently a plain list) with the following items:
  \item{distance }{The computed distance}
  \item{index1 }{Matched elements: indices in $x$}
  \item{index2 }{Corresponding mapped indices in $y$}
  \item{stepPatterns }{The step patterns used in the computation, in
    the matrix form documented above.}
  \item{costMatrix}{ If \code{keep.internals=TRUE}, the cumulative
    cost matrix}
  \item{directionMatrix}{ If \code{keep.internals=TRUE}, the
    directions of steps that would be taken at each alignment pair
    (integers indexing rows in stepPatterns)}
  }
\references{ TODO }
\author{Toni Giorgino }
\note{Cost matrices (both input and output) have query elements row-wise
  (first index), and template elements column-wise (second index). They
  print as usual, with indexes increasing down and right.  In some
  DTW  tutorials cost matrices follow  plot-like conventions
  (e.g.  one index grows upward). This may
  be confusing.
}
\seealso{  \code{\link[pkg:analogue]{distance}},
  \code{\link{dtwWindowingFunctions}}, \code{\link{outer}} }
\examples{
## A noisy sine wave as query
idx<-seq(0,6.28,len=100)
query<-sin(idx)+runif(100)/10;

## A cosine is for template; they are offset by 25 samples
template<-cos(idx)
plot(template); lines(query);

## Find the best match (approx 1sec on Pentium 4)
dtw(query,template)->alignment;

## Display the mapping
plot(alignment$index1,alignment$index2);

## That's all: 25 samples off-diagonal alignment
lines(1:100-25,col="red")


}

\concept{Dynamic Time Warp}
\concept{Align timeseries}
\concept{Minimum cumulative cost}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ ts }
\keyword{ optimize }


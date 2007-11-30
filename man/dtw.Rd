\name{dtw}
\alias{dtw}
\alias{is.dtw}



\title{Dynamic Time Warp}
\description{
  Compute Dynamic Time Warp
  and find optimal alignment between two time series.
}
\usage{
dtw(x, y=NULL, distance.function=euclideanSquared, step.pattern="s",
  window.type="none", keep.internals=FALSE, ...)

is.dtw(d)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ query vector OR local cost matrix }
  \item{y}{ template vector, or unused if cost matrix given }
  %  \item{partial}{ ~~Describe \code{partial} here~~ }
  \item{distance.function}{ pointwise distance function. See e.g. \code{\link{dtwDistanceFunctions}} }
  \item{step.pattern}{ step pattern. Char: "s" (symmetric), "a"
    (asymmetric), or step matrix containing the allowed
    steps with their cost (see \code{\link{stepPattern}})}
  \item{window.type}{ windowing function, character or function. Char: "none",
    "itakura", "sakoechiba", "slantedband". Function: boolean of two arguments.}
  \item{keep.internals}{don't discard the cumulative cost matrix and other
    internal structures}
  \item{d}{an arbitrary R object}
  \item{...}{additional arguments, passed to \code{window.function}}
}

\details{

  The function performs Dynamic Time Warp (DTW) and computes the optimal
  alignment between two time series \code{x} and \code{y}, given as
  numeric vectors.  The ``optimal'' alignment minimizes the sum of
  distances between aligned elements. Lengths of \code{x} and \code{y} may be different.
  
  The (local) distance between elements of \code{x} (query) and \code{y}
  (template) is computed through the supplied \code{distance.function}
  (see \code{\link{dtwDistanceFunctions}}). 

  Multivariate time series and complex distance functions can be handled
  by supplying a local-distance
  matrix. Element \code{[i,j]} of the
  local-distance matrix is understood as the distance between element  \code{x[i]} and
  \code{y[j]}. The distance matrix has therefore \code{n=length(x)} rows and
  \code{m=length(y)} columns (see note below). 
  To generate a local distance matrix, the functions \code{outer}
  and \code{\link[pkg:analogue]{distance}} in package
  \pkg{analogue} may come handy.

  Several common variants of DTW are supported via the \code{step.pattern}
  argument, which defaults to \code{symmetric1} (White-Neely). Most common
  step patterns are supported and pre-defined; the 
  user can also write their own. See 
  \code{\link{stepPattern}} for details.

  Windowing is supported by supplying a name into the \code{window.type}
  argument (abbreviations allowed)  between the built-in types:

  \itemize{
	\item{\code{"none"}}{No windowing (default)}
	\item{\code{"sakoechiba"}}{A band around main diagonal}
	\item{\code{"slantedband"}}{A band around slanted diagonal}
	\item{\code{"itakura"}}{So-called Itakura parallelogram}
  }

  \code{window.type} can also be user-defined windowing function. 
  See \code{\link{dtwWindowingFunctions}} for
  all available windowing functions, details on user-defined windowing, 
  and a discussion of the (mis)naming of the "Itakura" parallelogram 
  as a global constraint. 

  Some windowing functions may require parameters, such as the
  \code{window.size} argument.

  A native (fast, compiled) version of the function is normally available.
  If it is not, an interpreted equivalent will be used as 
  a fall-back, with a warning.

  \code{is.dtw} tests whether the argument is of class \code{dtw}.

}



\value{
  An object of class \code{dtw} with the following items:
  \item{distance}{the computed distance \emph{not normalized}. Normalization 
	depends on the chosen step pattern.}
  \item{index1}{matched elements: indices in \code{x}}
  \item{index2}{corresponding mapped indices in \code{y}}
  \item{stepPatterns}{the \code{stepPattern} object used for the computation}
  \item{costMatrix}{if \code{keep.internals=TRUE}, the cumulative
    cost matrix}
  \item{directionMatrix}{if \code{keep.internals=TRUE}, the
    directions of steps that would be taken at each alignment pair
    (integers indexing step patterns)}
  }
  \references{
  Sakoe, H.; Chiba, S., \emph{Dynamic programming algorithm optimization for spoken word recognition,}
 Acoustics, Speech, and Signal Processing [see also IEEE Transactions on Signal Processing], IEEE Transactions on , 
vol.26, no.1, pp. 43-49, Feb 1978 URL: \url{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055} \cr
  }
\author{Toni Giorgino }
\note{Cost matrices (both input and output) have query elements row-wise
  (first index), and template elements column-wise (second index). They
  print according to the usual convention, with indexes increasing down and right.  Many
  DTW  papers and tutorials show matrices according to plot-like conventions, i.e.
  one index growing upwards. This may   be confusing.
}
\seealso{
  \code{\link{outer}},\code{\link[pkg:analogue]{distance}};
  \code{\link{dtwWindowingFunctions}},  windowing and global constraints;
  \code{\link{stepPattern}}, step patterns and local constraints
  including slope;
  \code{\link{plot.dtw}}, the plot method for DTW objects.
}


\examples{
## A noisy sine wave as query
idx<-seq(0,6.28,len=100);
query<-sin(idx)+runif(100)/10;

## A cosine is for template; sin and cos are offset by 25 samples
template<-cos(idx)
plot(template); lines(query);

## Find the best match (approx 1sec on Pentium 4)
## keep = TRUE so we can make a density plot later on
alignment<-dtw(query,template,keep=TRUE);


## Display the mapping
plot(alignment$index1,alignment$index2);

## That's all: 25 samples off-diagonal alignment
lines(1:100-25,col="red")


## Contour plots of the global cost


## A profile of the cumulative distance matrix
## similar to: plot(alignment,type="density") or dtwPlotDensity(alignment)

contour(alignment$costMatrix,col=terrain.colors(100),x=1:100,y=1:100,
	xlab="Query (noisy sine)",ylab="Template (cosine)");

lines(alignment$index1,alignment$index2,col="red",lwd=2);


## More plots on dtw.plot!


#########
## An hand-checkable example

ldist<-matrix(1,nrow=6,ncol=6);  # Matrix of ones
ldist[2,]<-0; ldist[,5]<-0;      # Mark a clear path of zeroes
ldist[2,5]<-.01;		 # Forcely cut the corner

ds<-dtw(ldist);			 # DTW with user-supplied local cost matrix
da<-dtw(ldist,step="a");	 # Also compute the asymmetric 
plot(ds$index1,ds$index2,pch=3);    	 # Symmetric: alignment follows the marked path
points(da$index1,da$index2,col="red");    	 # Asymmetric: visiting 1 is required twice

ds$distance;
da$distance;




}

\concept{Dynamic Time Warp}
\concept{Dynamic programming}
\concept{Align timeseries}
\concept{Minimum cumulative cost}



\keyword{ ts }
\keyword{ optimize }


\name{plot.dtw}
\alias{plot.dtw}
\alias{dtwPlot}
\alias{dtwPlotAlignment}
\alias{dtwPlotDensity}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting of dynamic time warp results}
\description{
  Methods for plotting dynamic time warp alignment objects returned
  by \code{\link{dtw}}.
}
\usage{

\method{plot}{dtw}(x, type="alignment", ...)

# an alias for dtw.plot
dtwPlot(x, type="alignment", ...)

dtwPlotAlignment(d, xlab="Query index", ylab="Template index", ...)
dtwPlotDensity(d, normalize="no",
		  xlab="Query index", ylab="Template index",
                  ...)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,d}{ \code{dtw} object, usually result of call to \code{\link{dtw}}}
  \item{xlab}{ label for the query axis}
  \item{ylab}{label for the template axis}
  \item{type}{alignment plot style}
  \item{normalize}{show per-step average cost instead of cumulative cost}
  
  \item{...}{additional arguments, passed to plotting functions}
}

\details{

  \code{dtwPlot} displays  alignment contained in \code{dtw} objects.

  Various plotting styles are available, passing strings to the
  \code{type} argument (may be abbreviated):

  \itemize{
    \item{\code{alignment}}{simple plot of the warping path}
    %    \item{\code{twoway}}{Point-by-point comparison}
    \item{\code{threeway}}{vis-a-vis visual inspection of the
       original timeseries and their alignment }
    \item{\code{density}}{show the cumulative cost matrix with the
      warping path overimposed}
  }

  Three-way plotting is documented separately in
  function \code{\link{dtwPlotThreeWay}}.

  \code{normalize} can be one of \code{"N"} of \code{"N+M"}.
  If set, \emph{average} cost per step is plotted instead
  of the cumulative one. Step averaging depends on the
  step pattern chosen. \code{N}  is suitable
  for asymmetric patterns and divides the distance
  by the number of steps along the query; \code{N+M}
  is for (normalizable) symmetric patterns,
  and divides by the Manhattan distance from
  the origin.
  
  Additional  parameters are carried on to the plotting
  functions: use with care.

}

\note{
  The density plot is more colorful than useful.
  }

\author{Toni Giorgino }

\seealso{
  \code{\link{dtwPlotThreeWay}} for details on  three-way plotting function.
}


\examples{
## Same example as in dtw

idx<-seq(0,6.28,len=100);
query<-sin(idx)+runif(100)/10;
template<-cos(idx)

alignment<-dtw(query,template,keep=TRUE);


## A profile of the cumulative distance matrix
## Contour plot of the global cost

dtwPlotDensity(alignment, main="Sine/cosine: symmetric  alignment, no constraints")



######

## A study of the "itakura" parallelogram


## A widely held misconception is that the "Itakura parallelogram" 
## (as described in the original article) is a global constraint.
## Instead, it arises from local slope restrictions. Anyway, an "itakuraWindow",
## is provided in this package. A comparison between the two follows.


## The local constraint: three sides of the parallelogram are seen
dtw(query,template,keep=TRUE,step=asymmetricItakura)->ita;
dtwPlot(ita,type="density",main="Slope-limited asymmetric step (Itakura)")

## Symmetric step with global parallelogram-shaped constraint
## Note how long (>2 steps) horizontal stretches are allowed within the window.
dtw(query,template,keep=TRUE,window=itakuraWindow)->ita;
dtwPlot(ita,type="density",main="Symmetric step with Itakura parallelogram window")


}

\concept{Dynamic Time Warp}


\keyword{ ts }
\keyword{ hplot }

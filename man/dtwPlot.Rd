\name{dtwPlot}
\alias{dtwPlot}
\alias{dtwPlotAlignment}
\alias{dtwPlotDensity}
\alias{plot.dtw}



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

dtwPlotAlignment(d, xlab="Query index", ylab="Reference index", ...)
dtwPlotDensity(d, normalize=FALSE,
		  xlab="Query index", ylab="Reference index",
                  ...)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,d}{ \code{dtw} object, usually result of call to \code{\link{dtw}}}
  \item{xlab}{ label for the query axis}
  \item{ylab}{label for the reference axis}
  \item{type}{alignment plot style}
  \item{normalize}{show per-step average cost instead of cumulative cost}
  
  \item{...}{additional arguments, passed to plotting functions}
}

\details{

  \code{dtwPlot} displays  alignment contained in \code{dtw} objects.

  Various plotting styles are available, passing strings to the
  \code{type} argument (may be abbreviated):

  \itemize{
    \item{\code{alignment}}{plots the warping curve in \code{d}}
    \item{\code{twoway}}{plots a point-by-point comparison, with matching lines}
    \item{\code{threeway}}{vis-a-vis inspection of the
        timeseries and their warping curve }
    \item{\code{density}}{displays the cumulative cost landscape with the
      warping path overimposed}
  }

  
  For two-way plotting, see documentation for function
  \code{\link{dtwPlotTwoWay}}.

  For three-way plotting, see documentation for function
  \code{\link{dtwPlotThreeWay}}.


  If \code{normalize} is \code{TRUE}, the \emph{average} cost per step
  is plotted instead of the cumulative one. Step averaging depends on
  the \code{\link{stepPattern}} used.
  
  Additional  parameters are carried on to the plotting
  functions: use with care.

}


\note{
  The density plot is more colorful than useful.
  }

\section{Warning}{ These functions are incompatible with mechanisms for
     arranging plots on a device: \code{par(mfrow)}, \code{layout} and
     \code{split.screen}.}

\author{Toni Giorgino }

\seealso{
  \code{\link{dtwPlotTwoWay}} for details on  two-way plotting function.
  \code{\link{dtwPlotThreeWay}} for details on  three-way plotting function.
}


\examples{
## Same example as in dtw

idx<-seq(0,6.28,len=100);
query<-sin(idx)+runif(100)/10;
reference<-cos(idx)

alignment<-dtw(query,reference,keep=TRUE);


## A profile of the cumulative distance matrix
## Contour plot of the global cost

dtwPlotDensity(alignment,
  main="Sine/cosine: symmetric  alignment, no constraints")



######
##
## A study of the "Itakura" parallelogram
##
## A widely held misconception is that the "Itakura parallelogram" (as
## described in the original article) is a global constraint.  Instead,
## it arises from local slope restrictions. Anyway, an "itakuraWindow",
## is provided in this package. A comparison between the two follows.


## The local constraint: three sides of the parallelogram are seen

dtw(query,reference,keep=TRUE,step=typeIIIc)->ita;
dtwPlot(ita,type="density",
        main="Slope-limited asymmetric step (Itakura)")

## Symmetric step with global parallelogram-shaped constraint. Note how
## long (>2 steps) horizontal stretches are allowed within the window.

dtw(query,reference,keep=TRUE,window=itakuraWindow)->ita;
dtwPlot(ita,type="density",
        main="Symmetric step with Itakura parallelogram window")


}

\concept{Dynamic Time Warp}
\concept{Warping function}


\keyword{ ts }
\keyword{ hplot }

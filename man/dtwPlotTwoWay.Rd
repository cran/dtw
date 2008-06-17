\name{dtwPlotTwoWay}
\alias{dtwPlotTwoWay}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting of dynamic time warp results: pointwise comparison}
\description{
  Display the query and template time series and their alignment,
  arranged for visual inspection.
}

\usage{
dtwPlotTwoWay(d,xts=NULL,yts=NULL, offset=0,
			ts.type="l",pch=21, 
                        match.indices=NULL,
			match.col="gray70",
			xlab="Index", ylab="Query value", 
			... )
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ an alignment result, object of class \code{dtw}}
  \item{xts}{ query vector }
  \item{yts}{  template vector }
  \item{xlab,ylab}{axis labels}
  \item{offset}{displacement between the timeseries, summed to template}
  \item{match.col}{color of the match lines}
  \item{match.indices}{indices for which to draw a visual
    guide}
  \item{ts.type,pch}{graphical parameters for timeseries plotting, passed to \code{matplot}}
  \item{...}{additional arguments, passed to \code{matplot}}
}

\details{

  The two vectors are displayed via the \code{\link{matplot}} functions;
  their appearance can be customized via the \code{type} and \code{pch}
  arguments (constants or vectors of two elements).  If \code{offset} is
  set, the template is shifted vertically by the given amount; this will
  be reflected by the \emph{right-hand} axis.

  If a vector of integers \code{match.indices} is provided, a visual
  guide is only drawn for the corresponding indices in the warping curve
  (match lines).

  If \code{xts} and \code{yts} are not supplied, they will be recovered
  from \code{d}, as long as it was created with the two-argument call of
  \code{\link{dtw}}  with \code{keep.internals=T}.  Only
  single-variate time series can be plotted this way.
  
}


\author{Toni Giorgino }

\note{
  When \code{offset} is set values on the left axis only apply to
  the query.
}

\section{Warning}{ The function is incompatible with mechanisms for
     arranging plots on a device: \code{par(mfrow)}, \code{layout} and
     \code{split.screen}.}


\seealso{ \code{\link{dtwPlot}} for other dtw plotting functions,
\code{\link{matplot}} for graphical parameters. }

\examples{

## A noisy sine wave as query
## A cosine is for template; sin and cos are offset by 25 samples

idx<-seq(0,6.28,len=100);
query<-sin(idx)+runif(100)/10;
template<-cos(idx)
dtw(query,template,step=asymmetricP1)->alignment;


## Equivalent to plot(alignment,xts=query,yts=template,type="two");
dtwPlotTwoWay(alignment,xts=query,yts=template);

## Beware of the template's y axis, may be confusing
plot(alignment,xts=query,yts=template,offset=-2,type="two");

legend("topright",c("Query","Template"), pch=21, col=1:6)



}

\concept{Dynamic Time Warp}


\keyword{ hplot }


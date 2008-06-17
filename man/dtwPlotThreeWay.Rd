\name{dtwPlotThreeWay}
\alias{dtwPlotThreeWay}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting of dynamic time warp results: annotated warping function}
\description{
  Display the query and template time series and their warping curve,
  arranged for visual inspection. 
}

\usage{

dtwPlotThreeWay(d,xts=NULL,yts=NULL,
                type.align="l",type.ts="l",
                match.indices=NULL,
	        margin=4, inner.margin=0.2, title.margin=1.5,
	        xlab="Query index",ylab="Template index",main="Timeseries alignment",
			... )

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ an alignment result, object of class \code{dtw}}
  \item{xts}{ query vector }
  \item{yts}{  template vector }
  \item{xlab}{ label for the query axis}
  \item{ylab}{label for the template axis}
  \item{main}{main title}
  \item{type.align}{line style for warping curve plot}
  \item{type.ts}{line style for timeseries plot}
  \item{match.indices}{indices for which to draw a visual
    guide}
  \item{margin}{outer figure margin }
  \item{inner.margin}{inner figure margin }
  \item{title.margin}{space on the top of figure}
  
  \item{...}{additional arguments, used for the warping curve}
}

\details{

  The query time series is plotted in the bottom panel, with indices
  growing rightwards and values upwards. Template is in the left panel,
  indices growing upwards and values leftwards. The warping curve panel
  matches indices, and therefore element (1,1) will be at the lower
  left, (N,M) at the upper right.

  If a vector of integers \code{match.indices} is provided, a visual
  guide is drawn for the corresponding indices in the warping curve
  (match lines).

  If \code{xts} and \code{yts} are not supplied, they will be recovered
  from \code{d}, as long as it was created with the two-argument call of
  \code{\link{dtw}} with \code{keep.internals=T}.  Only single-variate
  time series can be plotted.

  
}



\author{Toni Giorgino }


\section{Warning}{ The function is incompatible with mechanisms for
     arranging plots on a device: \code{par(mfrow)}, \code{layout} and
     \code{split.screen}. Appearance of the match lines and timeseries
     currently  can not be customized. }

\examples{

## A noisy sine wave as query
## A cosine is for template; sin and cos are offset by 25 samples

idx<-seq(0,6.28,len=100);
query<-sin(idx)+runif(100)/10;
template<-cos(idx)
dtw(query,template,keep=TRUE)->alignment;


## Beware of the template's y axis, may be confusing
## Equivalent to plot(alignment,xts=query,yts=template,type="three");
dtwPlotThreeWay(alignment,xts=query,yts=template);




}

\concept{Dynamic Time Warp}


\keyword{ hplot }


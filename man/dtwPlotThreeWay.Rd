\name{dtwPlotThreeWay}
\alias{dtwPlotThreeWay}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting of dynamic time warp results: annotated warping function}
\description{
  Plot a DTW results for visual inspection.
}

\usage{

dtwPlotThreeWay(d,xts=NULL,yts=NULL,type.align="p",type.ts="l",
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
  \item{type.align}{line style for alignment path plot}
  \item{type.ts}{line style for timeseries plot}
  \item{margin}{outer figure margin }
  \item{inner.margin}{inner figure margin }
  \item{title.margin}{space on the top of figure}
  
  \item{...}{additional arguments, passed to plotting functions}
}

\details{

  Display the query and template time series and their alignment,
  arranged so that alignment can be visually inspected.

  Only single-variate time series are supported, and they must be passed
  to the \code{xts} and \code{yts} arguments (query and template are not
  stored in the \code{dtw} object).
  
}


\author{Toni Giorgino }

\note{
  The value axis for the template time series grows rightwards.
  This may be counter-intuitive.
}

\section{Warning}{ These functions are incompatible with mechanisms for
     arranging plots on a device: \code{par(mfrow)}, \code{layout} and
     \code{split.screen}.}

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


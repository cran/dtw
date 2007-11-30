\name{dtwWindowingFunctions}
\alias{noWindow}
\alias{sakoeChibaWindow}
\alias{itakuraWindow}
\alias{dtwWindowingFunctions}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sakoe-Chiba band and other windowing for DTW}
\description{
  Windowing functions used by \code{\link{dtw}}.
}
\usage{
noWindow(iii, jjj);
sakoeChibaWindow(iii, jjj);
itakuraWindow(iii, jjj);
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iii}{ Index in the query }
  \item{jjj}{ Index in the template }
}
\details{

  These functions, and user-defined equivalents to be passed to the
  \code{window.type} argument in \code{\link{dtw}}, take two integer
  arguments and return a boolean value stating if the coordinates fall
  within the allowed window for warping paths.
  
  User-defined functions can read variables \code{template.size},
  \code{query.size} and \code{window.size}; these are forwarded by
  \code{dtw}. User-defined functions are free to implement any window
  shape, as long as at least one path is allowed between the initial and
  final alignment points, i.e., they are compatible with the DTW
  constraints.

  Windowing in the rectangular case is supported as a package-specific
  extension (see note).  If the cost matrix is non-square, the
  \emph{built-in} window types \code{sakoeChibaWindow} and
  \code{itakuraWindow} are automatically adjusted (see note). The
  Sakoe-Chiba band, in particular, is centered about the (jagged) line
  which joins element \code{[1,1]} to element \code{[n,m]}, and will be
  \code{window.size} columns wide.


  
}
\value{
  The functions return \code{TRUE}  if the given coordinates fall within the chosen
  warping window,  \code{FALSE} otherwise. User-defined functions should
  do the same.
}
\references{ TODO }
\author{Toni Giorgino}

\note{Allowing windowing in the rectangular case (\code{template.size !=
  query.size}) is a package-specific generalization. The Sakoe-Chiba
  band width is measured along the template (column-wise). }

\examples{

## The sakoeChibaWindow function is currently crudely defined as
function(iii,jjj) {
  diagj<-(iii*template.size/query.size);
  return(abs(jjj-diagj)<=window.size);
  }
}


\concept{Dynamic Time Warp}
\concept{Windowing}
\concept{Sakoe-Chiba Band}
\concept{Itakura Parallelogram}


\keyword{ ts }
    
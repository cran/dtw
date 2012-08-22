\name{dtw-package}
\alias{dtw-package}

\docType{package}
\title{
Dynamic Time Warp algorithms in R
}
\description{
  Dynamic Time Warp: find the optimal alignment between two time series.
}

\details{
\tabular{ll}{
Package: \tab  dtw		\cr
Type: 	 \tab  Package		\cr
Version: \tab  1.15		\cr
Date: 	 \tab  2012-8-22	\cr
License: \tab  GPL-2		\cr
}

Comprehensive implementation of Dynamic Time Warping (DTW) algorithms in
R.

The basic DTW algorithm computes the time axis stretch which optimally
maps one timeseries (query) onto another (reference); it outputs the
remaining cumulative distance between the two. DTW is widely used
e.g. for classification and clustering tasks in econometrics,
chemometrics and general timeseries mining.

The R implementation in dtw provides: 

\itemize{
\item    arbitrary windowing functions (global constraints), eg. the Sakoe-Chiba band; see \code{\link{dtwWindowingFunctions}}
\item    arbitrary transition types (also known as step patterns, slope constraints, local constraints, or DP-recursion rules). This includes dozens of well-known types; see \code{\link{stepPattern}}:
\itemize{
 \item        all step patterns classified by Rabiner-Juang, Sakoe-Chiba, and Rabiner-Myers;
 \item        symmetric and asymmetric;
 \item        Rabiner's smoothed variants;
 \item        arbitrary, user-defined slope constraints
}
\item    partial matches: open-begin, open-end, substring matches 
\item    proper, pattern-dependent, normalization (exact average distance per step)
\item    the Minimum Variance Matching (MVM) algorithm (Latecki et al.)
}

Multivariate timeseries can be aligned with arbitrary local distance
definitions, leveraging the \code{\link[proxy]{dist}} function of
package \pkg{proxy}. DTW itself becomes a distance function with the
dist semantics.

In addition to computing alignments, the package provides:
\itemize{
\item    methods for plotting alignments and warping functions in several classic styles (see plot gallery);
\item    graphical representation of step patterns;
\item    functions for applying a warping function, either direct or inverse;
\item    both fast native (C) and interpreted (R) cores.
}


Please see documentation for function \code{\link{dtw}}, which is the
main entry point to the package.

If you use this software, please cite it according to
\code{citation("dtw")}.  The package home page is at
\url{http://dtw.r-forge.r-project.org}.

To get the latest stable version from CRAN, use
\code{install.packages("dtw")}.  To get the development version
(possibly unstable), use
\code{install.packages("dtw",repos="http://r-forge.r-project.org")}.

}
\author{
  Toni Giorgino, Copyright (c) 2007-2012
  Istituto di Ingegneria Biomedica (ISIB-CNR)
  National Research Council of Italy

  Maintainer:  toni.giorgino@gmail.com
}

\seealso{ \code{\link{dtw}} for the main entry point to the package;
  \code{\link{dtwWindowingFunctions}} for global constraints;
  \code{\link{stepPattern}} for local constraints;
  \code{\link[analogue]{distance}},  \code{\link{outer}} for
  building a local cost matrix with multivariate
  timeseries and custom distance functions. }

\references{
  Toni Giorgino. \emph{Computing and Visualizing Dynamic Time Warping
        Alignments in R: The dtw Package.}  Journal of Statistical
        Software, 31(7), 1-24. \url{http://www.jstatsoft.org/v31/i07/}

  Rabiner, L. R., & Juang, B.-H. (1993). Chapter 4 in
  \emph{Fundamentals of speech
  recognition.} Englewood Cliffs, NJ: Prentice Hall.
}

\examples{
 library(dtw);
 ## demo(dtw);
}


\keyword{ package }
\keyword{ ts }

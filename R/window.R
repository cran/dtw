###############################################################
#                                                             #
#   Author: Toni Giorgino <toni.giorgino,gmail.com>           #
#       Istituto di Neuroscienze (IN-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id$
#                                                             #
###############################################################

##
## Document in dtwWindowingFunctions.Rd
##




#' Global constraints and windowing functions for DTW
#' 
#' Various global constraints (windows) which can be applied to the
#' \code{window.type} argument of \code{\link{dtw}}, including the Sakoe-Chiba
#' band, the Itakura parallelogram, and custom functions.
#' 
#' 
#' Windowing functions can be passed to the \code{window.type} argument in
#' \code{\link{dtw}} to put a global constraint to the warping paths allowed.
#' They take two integer arguments (plus optional parameters) and must return a
#' boolean value \code{TRUE} if the coordinates fall within the allowed region
#' for warping paths, \code{FALSE} otherwise.
#' 
#' User-defined functions can read variables \code{reference.size},
#' \code{query.size} and \code{window.size}; these are pre-set upon invocation.
#' Some functions require additional parameters which must be set (e.g.
#' \code{window.size}).  User-defined functions are free to implement any
#' window shape, as long as at least one path is allowed between the initial
#' and final alignment points, i.e., they are compatible with the DTW
#' constraints.
#' 
#' The \code{sakoeChibaWindow} function implements the Sakoe-Chiba band, i.e.
#' \code{window.size} elements around the \code{main} diagonal. If the window
#' size is too small, i.e. if \code{reference.size}-\code{query.size} >
#' \code{window.size}, warping becomes impossible.
#' 
#' An \code{itakuraWindow} global constraint is still provided with this
#' package.  See example below for a demonstration of the difference between a
#' local the two.
#' 
#' The \code{slantedBandWindow} (package-specific) is a band centered around
#' the (jagged) line segment which joins element \code{[1,1]} to element
#' \code{[query.size,reference.size]}, and will be \code{window.size} columns
#' wide. In other words, the "diagonal" goes from one corner to the other of
#' the possibly rectangular cost matrix, therefore having a slope of
#' \code{M/N}, not 1.
#' 
#' \code{dtwWindow.plot} visualizes a windowing function. By default it plots a
#' 200 x 220 rectangular region, which can be changed via \code{reference.size}
#' and \code{query.size} arguments.
#' 
#' @name dtwWindowingFunctions
#' @aliases noWindow sakoeChibaWindow slantedBandWindow itakuraWindow
#' dtwWindowingFunctions dtwWindow.plot
#' @export noWindow sakoeChibaWindow slantedBandWindow itakuraWindow dtwWindow.plot
#' @param iw index in the query (row) -- automatically set
#' @param jw index in the reference (column) -- automatically set
#' @param query.size size of the query time series -- automatically set
#' @param reference.size size of the reference time series -- automatically set
#' @param window.size window size, used by some windowing functions -- must be
#' set
#' @param fun a windowing function
#' @param ... additional arguments passed to windowing functions
#' @return Windowing functions return \code{TRUE} if the coordinates passed as
#' arguments fall within the chosen warping window, \code{FALSE} otherwise.
#' User-defined functions should do the same.
#' @note Although \code{dtwWindow.plot} resembles object-oriented notation,
#' there is not a such a dtwWindow class currently.
#' 
#' A widely held misconception is that the "Itakura parallelogram" (as
#' described in reference [2]) is a \emph{global} constraint, i.e. a window.
#' To the author's knowledge, it instead arises from the local slope
#' restrictions imposed to the warping path, such as the one implemented by the
#' \code{\link{typeIIIc}} step pattern.
#' @author Toni Giorgino
#' @references [1] Sakoe, H.; Chiba, S., \emph{Dynamic programming algorithm
#' optimization for spoken word recognition,} Acoustics, Speech, and Signal
#' Processing [see also IEEE Transactions on Signal Processing], IEEE
#' Transactions on , vol.26, no.1, pp. 43-49, Feb 1978 URL:
#' \url{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055} \cr \cr
#' [2] Itakura, F., \emph{Minimum prediction residual principle applied to
#' speech recognition,} Acoustics, Speech, and Signal Processing [see also IEEE
#' Transactions on Signal Processing], IEEE Transactions on , vol.23, no.1, pp.
#' 67-72, Feb 1975. URL:
#' \url{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1162641}
#' @keywords ts
#' @examples
#' 
#' 
#' ## Display some windowing functions
#' dtwWindow.plot(itakuraWindow, main="So-called Itakura parallelogram window")
#' dtwWindow.plot(slantedBandWindow, window.size=2,
#'   reference=13, query=17, main="The slantedBandWindow at window.size=2")
#' 
#' 
#' ## Asymmetric step with Sakoe-Chiba band
#' 
#' idx<-seq(0,6.28,len=100); 
#' query<-sin(idx)+runif(100)/10;
#' reference<-cos(idx);
#' 
#' asyband<-dtw(query,reference,keep=TRUE,
#'              step=asymmetric,
#'              window.type=sakoeChibaWindow,
#'              window.size=30                  );
#' 
#' dtwPlot(asyband,type="density",main="Sine/cosine: asymmetric step, S-C window")
#' 
#' 
#' 
NULL

## no warping window: no restrictions

`noWindow` <-
function(iw,jw,...) {
  return(TRUE);
}



## A band around the diagonal. The band includes the diagonal +-
## window.size.

`sakoeChibaWindow` <-
function(iw,jw,window.size,...) {
  return(abs(jw-iw)<=window.size);
}


## A band around the diagonal. The band includes the segment
## connecting [1,1] to [query.size,reference.size] window.size,
## measured along the second index (columns)

`slantedBandWindow` <-
function(iw,jw,query.size,reference.size,window.size,...) {
  diagj<-(iw*reference.size/query.size);
  return(abs(jw-diagj)<=window.size);
}



## "Itakura" parallelogram: see documentation.
## 

`itakuraWindow` <- 
function(iw,jw,query.size,reference.size,...) {
	## Shorthands
  	n<-query.size; 	
  	m<-reference.size;

	ok<- 	(jw <  2*iw) &
 		(iw <= 2*jw) &
		(iw >= n-1-2*(m-jw)) &
		(jw >  m-1-2*(n-iw)) ;

	return(ok);
}


## Plot a sample of the windowing function

dtwWindow.plot <- function(fun,query.size=200,reference.size=220,...) {
	n<-query.size;
  	m<-reference.size;

	mm<-matrix(0,n,m);
	mm[fun(row(mm),col(mm),query.size=n,reference.size=m,...)]<-1;

	image(  z=mm,
		x=1:n,
		y=1:m,
		xlab=sprintf("Query: samples 1..%d",n),
		ylab=sprintf("Reference: samples 1..%d",m)
	 );
}


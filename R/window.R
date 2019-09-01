
##
## Copyright (c) 2006-2019 of Toni Giorgino
##
## This file is part of the DTW package.
##
## DTW is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## DTW is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
## or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
## License for more details.
##
## You should have received a copy of the GNU General Public License
## along with DTW.  If not, see <http://www.gnu.org/licenses/>.
##

##
## Document in dtwWindowingFunctions.Rd
##

#' Global constraints and windowing functions for DTW
#' 
#' Various global constraints (windows) which can be applied to the
#' `window.type` argument of [dtw()], including the Sakoe-Chiba
#' band, the Itakura parallelogram, and custom functions.
#' 
#' 
#' Windowing functions can be passed to the `window.type` argument in
#' [dtw()] to put a global constraint to the warping paths allowed.
#' They take two integer arguments (plus optional parameters) and must return a
#' boolean value `TRUE` if the coordinates fall within the allowed region
#' for warping paths, `FALSE` otherwise.
#' 
#' User-defined functions can read variables `reference.size`,
#' `query.size` and `window.size`; these are pre-set upon invocation.
#' Some functions require additional parameters which must be set (e.g.
#' `window.size`).  User-defined functions are free to implement any
#' window shape, as long as at least one path is allowed between the initial
#' and final alignment points, i.e., they are compatible with the DTW
#' constraints.
#' 
#' The `sakoeChibaWindow` function implements the Sakoe-Chiba band, i.e.
#' `window.size` elements around the `main` diagonal. If the window
#' size is too small, i.e. if `reference.size`-`query.size` >
#' `window.size`, warping becomes impossible.
#' 
#' An `itakuraWindow` global constraint is still provided with this
#' package.  See example below for a demonstration of the difference between a
#' local the two.
#' 
#' The `slantedBandWindow` (package-specific) is a band centered around
#' the (jagged) line segment which joins element `[1,1]` to element
#' `[query.size,reference.size]`, and will be `window.size` columns
#' wide. In other words, the "diagonal" goes from one corner to the other of
#' the possibly rectangular cost matrix, therefore having a slope of
#' `M/N`, not 1.
#' 
#' `dtwWindow.plot` visualizes a windowing function. By default it plots a
#' 200 x 220 rectangular region, which can be changed via `reference.size`
#' and `query.size` arguments.
#' 
#' @name dtwWindowingFunctions
#' @aliases noWindow sakoeChibaWindow slantedBandWindow itakuraWindow
#' dtwWindowingFunctions dtwWindow.plot
#' @param iw index in the query (row) -- automatically set
#' @param jw index in the reference (column) -- automatically set
#' @param query.size size of the query time series -- automatically set
#' @param reference.size size of the reference time series -- automatically set
#' @param window.size window size, used by some windowing functions -- must be
#' set
#' @param fun a windowing function
#' @param ... additional arguments passed to windowing functions
#' @return Windowing functions return `TRUE` if the coordinates passed as
#' arguments fall within the chosen warping window, `FALSE` otherwise.
#' User-defined functions should do the same.
#' @note Although `dtwWindow.plot` resembles object-oriented notation,
#' there is not a such a dtwWindow class currently.
#' 
#' A widely held misconception is that the "Itakura parallelogram" (as
#' described in reference 2) is a *global* constraint, i.e. a window.
#' To the author's knowledge, it instead arises from the local slope
#' restrictions imposed to the warping path, such as the one implemented by the
#' [typeIIIc()] step pattern.
#' @author Toni Giorgino
#' @references 
#' 1. Sakoe, H.; Chiba, S., *Dynamic programming algorithm
#' optimization for spoken word recognition,* Acoustics, Speech, and Signal
#' Processing, IEEE
#' Transactions on , vol.26, no.1, pp. 43-49, Feb 1978 URL:
#' <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055> 
#' 2. Itakura, F., *Minimum prediction residual principle applied to
#' speech recognition,* Acoustics, Speech, and Signal Processing, 
#' IEEE Transactions on , vol.23, no.1, pp.
#' 67-72, Feb 1975. URL:
#' <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1162641>
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




## A band around the diagonal. The band includes the diagonal +-
## window.size.
#' @export
#' @rdname dtwWindowingFunctions
`sakoeChibaWindow` <-
function(iw,jw,window.size,...) {
  return(abs(jw-iw)<=window.size);
}


## A band around the diagonal. The band includes the segment
## connecting [1,1] to [query.size,reference.size] window.size,
## measured along the second index (columns)
#' @export
#' @rdname dtwWindowingFunctions
`slantedBandWindow` <-
function(iw,jw,query.size,reference.size,window.size,...) {
  diagj<-(iw*reference.size/query.size);
  return(abs(jw-diagj)<=window.size);
}



## "Itakura" parallelogram: see documentation.
## 
#' @export
#' @rdname dtwWindowingFunctions
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
#' @export
#' @rdname dtwWindowingFunctions
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

## no warping window: no restrictions
#' @export
#' @rdname dtwWindowingFunctions
`noWindow` <-
    function(iw,jw,...) {
        return(TRUE);
    }



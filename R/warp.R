
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


## Apply the warping curve to a given timeseries. If the curve is
## multi-valued, average multiple points, with a warning. Argument
## "inverse=T" warps the inverse of the curve (ie a reference into a
## query).

## Direct:
##  all points in query have at least one image (NO?)
##  only one image, if asymmetric
##  if partial, image <= reference
##  x should be as long as the range in ix
##  there could be gaps in ix, to be interpolated,
##  depending on the step pattern

## for each point in the reference space, do an interpolated lookup
## into the y->x mapping

## Inverse: considering reference as domain
##  if partial, some trailing points may have no image
##  otherwise, all points have one or more images


## sortedXyData





#' Apply a warping to a given timeseries
#' 
#' 
#' Returns the indexing required to apply the optimal warping curve to a given
#' timeseries (warps either into a query or into a reference).
#' 
#' 
#' The warping is returned as a set of indices, which can be used to subscript
#' the timeseries to be warped (or rows in a matrix, if one wants to warp a
#' multivariate time series).  In other words, `warp` converts the warping
#' curve, or its inverse, into a function in the explicit form.
#' 
#' Multiple indices that would be mapped to a single point are averaged, with a
#' warning. Gaps in the index sequence are filled by linear interpolation.
#' 
#' @param d `dtw` object specifying the warping curve to apply
#' @param index.reference `TRUE` to warp a reference, `FALSE` to warp
#' a query
#' @return A list of indices to subscript the timeseries.
#' @author Toni Giorgino
#' @seealso Examples in [dtw()] show how to *graphically* apply
#' the warping via parametric plots.
#' @keywords ts
#' @examples
#' 
#' idx<-seq(0,6.28,len=100);
#' query<-sin(idx)+runif(100)/10;
#' reference<-cos(idx)
#' 
#' alignment<-dtw(query,reference);
#' 
#' 
#' wq<-warp(alignment,index.reference=FALSE);
#' wt<-warp(alignment,index.reference=TRUE);
#' 
#' old.par <- par(no.readonly = TRUE);
#' par(mfrow=c(2,1));
#' 
#' plot(reference,main="Warping query");
#'   lines(query[wq],col="blue");
#' 
#' plot(query,type="l",col="blue",
#'   main="Warping reference");
#'   points(reference[wt]);
#' 
#' par(old.par);
#' 
#' 
#' ##############
#' ##
#' ## Asymmetric step makes it "natural" to warp
#' ## the reference, because every query index has
#' ## exactly one image (q->t is a function)
#' ##
#' 
#' alignment<-dtw(query,reference,step=asymmetric)
#' wt<-warp(alignment,index.reference=TRUE);
#' 
#' plot(query,type="l",col="blue",
#'   main="Warping reference, asymmetric step");
#'   points(reference[wt]);
#' 
#' 
#' 
#' 
#' @export warp
warp <- function(d,index.reference=FALSE) {

  if(!is.dtw(d))
    stop("dtw object required");
  

  if(!index.reference) {
    ## warp QUERY into reference space
    iset<-d$index1;
    jset<-d$index2;
  } else {
    ## warp REFERENCE into query
    iset<-d$index2;
    jset<-d$index1;
  }
  jmax<-max(jset);

  ## rebuild  index, interpolating holes
  ii<-stats::approx(x=jset,y=iset,1:jmax);
  return(ii$y);    
     
}


## > ss<-warp.dtw(al,query)
## Warning message:
## In approx(x = jset, y = iset, 1:jmax) : collapsing to unique 'x' values
## > plot(reference);lines(query)
## > plot(reference);lines(ss)
## > st<-warp.dtw(al,reference,inverse=T)
## Warning message:
## In approx(x = jset, y = iset, 1:jmax) : collapsing to unique 'x' values
## > Cairo()
## > plot(query);lines(st)
## > plot(query,type="l");points(st)
## > 


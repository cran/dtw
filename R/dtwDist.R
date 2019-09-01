
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


## Compute a dissimilarity matrix, akin to "dist", analogue/distance,
## vegan/vegdist, etc. , based on the dtw "distance" measure.


## Apply FUN to all row pairs




#' Compute a dissimilarity matrix
#' 
#' Compute the dissimilarity matrix between a set of single-variate timeseries.
#' 
#' 
#' `dtwDist` computes a dissimilarity matrix, akin to [dist()],
#' based on the Dynamic Time Warping definition of a distance between
#' single-variate timeseries.
#' 
#' The `dtwDist` command is a synonym for the [proxy::dist()]
#' function of package \pkg{proxy}; the DTW distance is registered as
#' `method="DTW"` (see examples below).
#' 
#' The timeseries are stored as rows in the matrix argument `m`. In other
#' words, if `m` is an N * T matrix, `dtwDist` will build N*N ordered
#' pairs of timeseries, perform the corresponding N*N `dtw` alignments,
#' and return all of the results in a matrix. Each of the timeseries is T
#' elements long.
#' 
#' `dtwDist` returns a square matrix, whereas the `dist` object is
#' lower-triangular. This makes sense because in general the DTW "distance" is
#' not symmetric (see e.g.  asymmetric step patterns).  To make a square matrix
#' with the [proxy::dist()] function sematics, use the two-arguments
#' call as `dist(m,m)`. This will return a square `crossdist` object.
#' 
#' @param mx numeric matrix, containing timeseries as rows
#' @param my numeric matrix, containing timeseries as rows (for cross-distance)
#' @param ... arguments passed to the [dtw()] call
#' @return A square matrix whose element `[i,j]` holds the Dynamic Time
#' Warp distance between row `i` (query) and `j` (reference) of
#' `mx` and `my`, i.e.  `dtw(mx[i,],my[j,])$distance`.
#' @note To convert a square cross-distance matrix (`crossdist` object) to
#' a symmetric [dist()] object, use a suitable conversion strategy
#' (see examples).
#' @author Toni Giorgino
#' @seealso Other "distance" functions are: [dist()],
#' [vegan::vegdist()] in package `vegan`,
#' [analogue::distance()] in package `analogue`, etc.
#' @keywords ts
#' @examples
#' 
#' 
#' ## Symmetric step pattern => symmetric dissimilarity matrix;
#' ## no problem coercing it to a dist object:
#' 
#' m <- matrix(0,ncol=3,nrow=4)
#' m <- row(m)
#' dist(m,method="DTW");
#' 
#' # Old-fashioned call style would be:
#' #   dtwDist(m)
#' #   as.dist(dtwDist(m))
#' 
#' 
#' 
#' ## Find the optimal warping _and_ scale factor at the same time.
#' ## (There may be a better, analytic way)
#' 
#' # Prepare a query and a reference
#' 
#' query<-sin(seq(0,4*pi,len=100))
#' reference<-cos(seq(0,4*pi,len=100))
#' 
#' # Make a set of several references, scaled from 0 to 3 in .1 increments.
#' # Put them in a matrix, in rows
#' 
#' scaleSet <- seq(0.1,3,by=.1)
#' referenceSet<-outer(1/scaleSet,reference)
#' 
#' # The query has to be made into a 1-row matrix.
#' # Perform all of the alignments at once, and normalize the result.
#' 
#' dist(t(query),referenceSet,meth="DTW")->distanceSet
#' 
#' # The optimal scale for the reference is 1.0
#' plot(scaleSet,scaleSet*distanceSet,
#'   xlab="Reference scale factor (denominator)",
#'   ylab="DTW distance",type="o",
#'   main="Sine vs scaled cosine alignment, 0 to 4 pi")
#' 
#' 
#' 
#' 
#' 
#' ## Asymmetric step pattern: we can either disregard part of the pairs
#' ## (as.dist), or average with the transpose
#' 
#' mm <- matrix(runif(12),ncol=3)
#' dm <- dist(mm,mm,method="DTW",step=asymmetric); # a crossdist object
#' 
#' # Old-fashioned call style would be:
#' #   dm <- dtwDist(mm,step=asymmetric)
#' #   as.dist(dm)
#' 
#' 
#' ## Symmetrize by averaging:
#' (dm+t(dm))/2
#' 
#' 
#' ## check definition
#' stopifnot(dm[2,1]==dtw(mm[2,],mm[1,],step=asymmetric)$distance)
#' 
#' 
#' 
#' @export dtwDist
dtwDist <- function(mx,my=mx,...) {
  mye<-function(y,x,FUN,...) {
    apply(x,1,FUN,y,...);
  }

  apply(my,1,mye,mx,dtwpairdist,...);
}






###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino,gmail.com>           #
#       Istituto di Neuroscienze (IN-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id$
#                                                             #
###############################################################


## Compute a dissimilarity matrix, akin to "dist", analogue/distance,
## vegan/vegdist, etc. , based on the dtw "distance" measure.


## Apply FUN to all row pairs




#' Compute a dissimilarity matrix
#' 
#' Compute the dissimilarity matrix between a set of single-variate timeseries.
#' 
#' 
#' \code{dtwDist} computes a dissimilarity matrix, akin to \code{\link{dist}},
#' based on the Dynamic Time Warping definition of a distance between
#' single-variate timeseries.
#' 
#' The \code{dtwDist} command is a synonym for the \code{\link[proxy]{dist}}
#' function of package \pkg{proxy}; the DTW distance is registered as
#' \code{method="DTW"} (see examples below).
#' 
#' The timeseries are stored as rows in the matrix argument \code{m}. In other
#' words, if \code{m} is an N * T matrix, \code{dtwDist} will build N*N ordered
#' pairs of timeseries, perform the corresponding N*N \code{dtw} alignments,
#' and return all of the results in a matrix. Each of the timeseries is T
#' elements long.
#' 
#' \code{dtwDist} returns a square matrix, whereas the \code{dist} object is
#' lower-triangular. This makes sense because in general the DTW "distance" is
#' not symmetric (see e.g.  asymmetric step patterns).  To make a square matrix
#' with the \code{\link[proxy]{dist}} function sematics, use the two-arguments
#' call as \code{dist(m,m)}. This will return a square \code{crossdist} object.
#' 
#' @param mx numeric matrix, containing timeseries as rows
#' @param my numeric matrix, containing timeseries as rows (for cross-distance)
#' @param ... arguments passed to the \code{\link{dtw}} call
#' @return A square matrix whose element \code{[i,j]} holds the Dynamic Time
#' Warp distance between row \code{i} (query) and \code{j} (reference) of
#' \code{mx} and \code{my}, i.e.  \code{dtw(mx[i,],my[j,])$distance}.
#' @note To convert a square cross-distance matrix (\code{crossdist} object) to
#' a symmetric \code{\link{dist}} object, use a suitable conversion strategy
#' (see examples).
#' @author Toni Giorgino
#' @seealso Other "distance" functions are: \code{\link{dist}},
#' \code{\link[vegan]{vegdist}} in package \code{vegan},
#' \code{\link[analogue]{distance}} in package \code{analogue}, etc.
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





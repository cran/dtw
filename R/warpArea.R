
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


## Compute the (approximate) area between the warping function and the
## diagonal (in unit steps). 






#' Compute Warping Path Area
#' 
#' Compute the area between the warping function and the diagonal (no-warping)
#' path, in unit steps.
#' 
#' 
#' Above- and below- diagonal unit areas all count *plus* one (they do not
#' cancel with each other).  The "diagonal" goes from one corner to the other
#' of the possibly rectangular cost matrix, therefore having a slope of
#' `M/N`, not 1, as in [slantedBandWindow()].
#' 
#' The computation is approximate: points having multiple correspondences are
#' averaged, and points without a match are interpolated. Therefore, the area
#' can be fractionary.
#' 
#' @param d an object of class `dtw`
#' @return The area, not normalized by path length or else.
#' @note There could be alternative definitions to the area, including
#' considering the envelope of the path.
#' @author Toni Giorgino
#' @keywords ts
#' @examples
#' 
#'   ds<-dtw(1:4,1:8);
#' 
#'   plot(ds);lines(seq(1,8,len=4),col="red");
#' 
#'   warpArea(ds)
#' 
#'   ## Result: 6
#'   ##  index 2 is 2 while diag is 3.3  (+1.3)
#'   ##        3    3               5.7  (+2.7)
#'   ##        4   4:8 (avg to 6)    8   (+2  )
#'   ##                                 --------
#'   ##                                     6
#' 
#' @export warpArea
warpArea <- function(d) {

  if(!is.dtw(d))
    stop("dtw object required");
  
  ## rebuild query->templ map, interpolating holes
  ii<-stats::approx(x=d$index1,y=d$index2,1:d$N);

  dg<-seq(from=1,to=d$M,len=d$N);

  ad<-abs(ii$y-dg);
  sum(ad);
     
}




## Exmp:
##   t <- localWarpingStretch(alignment)
##
##   # local warping amount, based on the warping function
##   plot(t)
##
##   # lwa, based on the input position
##   plot(t~alignment$index1)
##
##   # or its pointwise approximation
##   plot(approx(alignment$index1,t,1:alignemnt$N)$y)

##
##   diff(localWarpingStretch) contains +1 for each "insertion" and -1
##    for each "deletion". Since we have max N deletions + M
##    insertions, a reasonable normalization would be to divide
##    sum(abs(diff())) by 2N
##
##   localWarpingCost is the amount of mismatch
##    at each point ("local substitution cost")



## Return how far from the diagonal is each point on the warping
## function. A good normalization factor can be 2 * N, so that maximum
## stretch is 1.

.localWarpingStretch <- function(d) {
    ## The diagonal line
    # dg <- seq(from = 1, to = d$M, len = d$N)

    ## get local copies
    id1 <- d$index1;
    id2 <- d$index2;

    ## remap reference indices to a square alignment
    id2 <- id2*d$N/d$M;

    ## return the local distance from the diagonal
    id1-id2;
}



## Return the local costs along the warping path i.e. d[i[t],j[t]] . A
## reasonable normalization could be to d$distance, so that each
## element would be the fraction of cost accumulated at that step.

.localWarpingCost <- function(d) {
    if(is.null(d$localCostMatrix))
        stop("A dtw object with keep.internals=TRUE is required");

    diag(d$localCostMatrix[d$index1,d$index2]);
}


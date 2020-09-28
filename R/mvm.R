
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


#' Minimum Variance Matching algorithm
#' 
#' Step patterns to compute the Minimum Variance Matching (MVM) correspondence
#' between time series
#' 
#' 
#' The Minimum Variance Matching algorithm (1) finds the non-contiguous parts
#' of reference which best match the query, allowing for arbitrarily long
#' "stretches" of reference to be excluded from the match. All elements of the
#' query have to be matched. First and last elements of the query are anchored
#' at the boundaries of the reference.
#' 
#' The `mvmStepPattern` function creates a `stepPattern` object which
#' implements this behavior, to be used with the usual [dtw()] call
#' (see example). MVM is computed as a special case of DTW, with a very large,
#' asymmetric-like step pattern.
#' 
#' The `elasticity` argument limits the maximum run length of reference
#' which can be skipped at once. If no limit is desired, set `elasticity`
#' to an integer at least as large as the reference (computation time grows
#' linearly).
#' 
#' @name mvm
#' @aliases mvm mvmStepPattern
#' @param elasticity integer: maximum consecutive reference elements skippable
#' @return A step pattern object.
#' @author Toni Giorgino
#' @seealso Other objects in [stepPattern()].
#' @references Latecki, L. J.; Megalooikonomou, V.; Wang, Q. & Yu, D.
#' *An elastic partial shape matching technique* Pattern Recognition,
#' 2007, 40, 3069-3080. \doi{10.1016/j.patcog.2007.03.004}
#' @keywords ts
#' @examples
#' 
#' 
#' ## The hand-checkable example given in Fig. 5, ref. [1] above
#' diffmx <- matrix( byrow=TRUE, nrow=5, c(
#'   0,  1,  8,  2,  2,  4,  8,
#'   1,  0,  7,  1,  1,  3,  7,
#'  -7, -6,  1, -5, -5, -3,  1,
#'  -5, -4,  3, -3, -3, -1,  3,
#'  -7, -6,  1, -5, -5, -3,  1 ) ) ;
#' 
#' ## Cost matrix
#' costmx <- diffmx^2;
#' 
#' ## Compute the alignment
#' al <- dtw(costmx,step.pattern=mvmStepPattern(10))
#' 
#' ## Elements 4,5 are skipped
#' print(al$index2)
#' 
#' plot(al,main="Minimum Variance Matching alignment")
#' 
#' 
#' 
#' @export
mvmStepPattern <- function(elasticity=20) {
  size <- elasticity;
  
  pn <- rep(1:size,each=2); # 112233...
  dx <- rep(c(1,0),size);   # 101010
  dy <- dx * ( (1:(2*size)) +1)/2;  # 102030
  w <-  rep(c(-1,1),size);
  
  tmp <- cbind(pn,dx,dy,w);
#  tmp <- as.vector(t(tmp));
  sp <- stepPattern(tmp);
  
  attr(sp,"norm") <- "N";
  attr(sp,"call") <- match.call();
  return(sp);
}





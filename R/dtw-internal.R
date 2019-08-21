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
## $Id$
##

## Internal functions for the dtw package.
## Not to be used by the user.



#' Internal dtw Functions
#' 
#' Internal dtw functions
#' 
#' These are not to be called by the user. Frontend to the DTW package is the
#' \code{\link{dtw}} function.
#' 
#' @name dtw-internal
#' @aliases backtrack dtwpairdist globalCostMatrix
#' @keywords internal
NULL

## Function applying dtw and only returning the
## distance
dtwpairdist <- function(...) {
  dtw(distance.only=TRUE,...)$distance;
}

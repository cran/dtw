
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
## $Id$
##

## Internal functions for the dtw package.
## Not to be used by the user.



#' Internal dtw Functions
#' 
#' Internal dtw functions
#' 
#' These are not to be called by the user. Frontend to the DTW package is the
#' [dtw()] function.
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


###############################################################
#                                                             #
#   Author: Toni Giorgino <toni.giorgino,gmail.com>           #
#       Istituto di Neuroscienze (IN-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id: dtw-internal.R 425 2016-08-25 19:48:58Z tonig $
#                                                             #
###############################################################

##
## $Id: dtw-internal.R 425 2016-08-25 19:48:58Z tonig $
##

## Internal functions for the dtw package.
## Not to be used by the user.


## Function applying dtw and only returning the
## distance
dtwpairdist <- function(...) {
  dtw(distance.only=TRUE,...)$distance;
}

###############################################################
#                                                             #
#   Author: Toni Giorgino <toni.giorgino@gmail.com>           #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: dtw-internal.R 94 2008-01-08 16:44:46Z tonig $
#                                                             #
###############################################################

##
## $Id: dtw-internal.R 94 2008-01-08 16:44:46Z tonig $
##

## Internal functions for the dtw package.
## Not to be used by the user.


## Function applying dtw and only returning the
## distance
dtwpairdist <- function(...) {
  dtw(distance.only=TRUE,...)$distance;
}

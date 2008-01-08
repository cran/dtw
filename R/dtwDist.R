###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino@gmail.com>           #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: distance.R 51 2007-12-11 10:59:08Z tonig $
#                                                             #
###############################################################


## Compute a dissimilarity matrix, akin to "dist", analogue/distance,
## vegan/vegdist, etc. , based on the dtw "distance" measure.


## Apply FUN to all row pairs
dtwDist <- function(m,...) {
  mye<-function(y,x,FUN,...) {
    apply(x,1,FUN,y,...);
  }

  apply(m,1,mye,m,dtwpairdist,...);
}





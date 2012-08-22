###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino,gmail.com>           #
#       Istituto di Ingegneria Biomedica (ISIB-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id: distance.R 51 2007-12-11 10:59:08Z tonig $
#                                                             #
###############################################################


## Compute a dissimilarity matrix, akin to "dist", analogue/distance,
## vegan/vegdist, etc. , based on the dtw "distance" measure.


## Apply FUN to all row pairs
dtwDist <- function(mx,my=mx,...) {
  mye<-function(y,x,FUN,...) {
    apply(x,1,FUN,y,...);
  }

  apply(my,1,mye,mx,dtwpairdist,...);
}





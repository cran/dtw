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


## Euclidean squared distance.
## should work on vectors (in R^n) too

## `euclideanSquared` <-
## function(a,b) {
##   z<-a-b;
##   z<-drop(z %*% z);                     #inner dot product
##   return (z);
## }

`euclideanSquared` <-
  function(a,b) {
    return((a-b)^2);
  }


`euclideanDistance` <-
  function(a,b) {
    return(abs(a-b));
  }

###############################################################
#                                                             #
#   Author: Toni Giorgino <toni.giorgino@gmail.com>           #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id$
#                                                             #
###############################################################




`dtw` <-
function(x, y=NULL,
         distance.function=euclideanSquared,
         step.pattern="s",
         window.type="none",
         window.size=10,
         keep.internals=FALSE ) {

  lm <- NULL;

  ## if matrix given
  if(is.matrix(x)) {
    lm <- x;
  } else {
    ## two ts. given TODO handle multivariate
    ## as.numeric handles ts, but not multivariates
    lm <- outer(as.numeric(x),
                as.numeric(y),
                FUN=distance.function);
  } 


  ## shorthand names
  n <- dim(lm)[1];
  m <- dim(lm)[2];

  gcm <- globalCostMatrix(lm, step.pattern, window.type, window.size);

  jmin <- m;

  ## result: distance (add to existing list gcm?)
  distance <- gcm$costMatrix[n,jmin];

  ## perform the backtrack
  mapping <- backtrack(jmin,gcm);

  ## append to existing list gcm, for now
  gcm$distance <- distance;
  gcm$index1 <- mapping$index1;
  gcm$index2 <- mapping$index2;

  ## delete sizey intermediate steps 
  if(!keep.internals) {
    gcm$costMatrix<-NULL;
    gcm$directionMatrix<-NULL;
  }

  ## if a dtw object is to be sponsored:
  # class(gcm) <- "dtw";

  return(gcm);
}


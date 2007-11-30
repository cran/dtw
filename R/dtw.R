###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino@gmail.com>           #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: dtw.R 69 2007-12-31 13:38:41Z tonig $
#                                                             #
###############################################################


##
## Frontend stuff, including coercing shorthand types
##


`dtw` <-
function(x, y=NULL,
         distance.function=euclideanSquared,
         step.pattern="s",
         window.type="none",
         keep.internals=FALSE,
         ... ) {

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


  ## Now we have a function
  wfun<-.canonicalizeWindowFunction(window.type);
  

  ## Now we have a step pattern
  dir<-.canonicalizeStepPattern(step.pattern);


  ## shorthand names
  n <- dim(lm)[1];
  m <- dim(lm)[2];

  if(is.loaded("computeCM")) {
    gcm <- globalCostNative(lm, step.matrix=dir,
                            window.function=wfun, ...);
  } else {
    warning("Native dtw implementation not available: using (slow) interpreted fallback");
    gcm <- globalCostMatrix(lm, step.matrix=dir,
                            window.function=wfun, ...);
  }


  
  jmin <- m;


  ## remember size
  gcm$N <- n;
  gcm$M <- m;
  
  ## result: distance (add to existing list gcm?)
  distance <- gcm$costMatrix[n,jmin];

  ## alignment valid?
  if(is.na(distance)) { stop("No warping paths exists that is allowed by costraints"); }

  ## perform the backtrack
  mapping <- backtrack(jmin,gcm);

  ## append to existing list gcm, for now
  ## perhaps replace by attr()
  gcm$distance <- distance;
  gcm$index1 <- mapping$index1;
  gcm$index2 <- mapping$index2;

  ## store call
  gcm$call <- match.call();

  ## delete sizey intermediate steps 
  if(!keep.internals) {
    gcm$costMatrix<-NULL;
    gcm$directionMatrix<-NULL;
  } 

  ## if a dtw object is to be sponsored:
  class(gcm) <- "dtw";

  return(gcm);
}


##############################
## OO class check
is.dtw <- function(d) {
    return(inherits(d,"dtw"));
}




## Replace  char window.type  with appropriate
## windowing FUNCTION 

.canonicalizeWindowFunction <- function(w) {
  if(is.function(w)) {
    return(w);
  }

  # else 
  wt<-pmatch(w,c("none","sakoechiba","itakura","slantedband"));
  if(is.na(wt)) {
    stop("Ambiguous or unsupported char argument for window.type");
  } 

  wfun<-switch(wt,
	noWindow,
	sakoeChibaWindow,
	itakuraWindow,
	slantedBandWindow);

  return(wfun);
}



## Replace char step.pattern with matrix equivalent
## Step patterns: \delta i, \delta j, cost
##                                      ...
## all deltas MUST be positive (otherwise we violate monotonicity)

.canonicalizeStepPattern <- function(step.pattern) {
  if(is.stepPattern(step.pattern))
    return(step.pattern);

  # else 
  sp<-pmatch(step.pattern,c("symmetric","asymmetric"));
  if(is.na(sp))
    stop("Ambiguous or unsupported char argument for step.pattern");

  dir<-switch(sp,
	symmetric1,
	asymmetric);

  return(dir);
}




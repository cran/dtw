###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino@gmail.com>           #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: dtw.R 83 2008-01-04 00:25:00Z tonig $
#                                                             #
###############################################################


##
## Frontend stuff, including coercing shorthand types
##


`dtw` <-
function(x, y=NULL,
         dist.method="Euclidean",
         step.pattern="s",
         window.type="none",
         keep.internals=FALSE,
         distance.only=FALSE,
         ... ) {

  lm <- NULL;


  
  ## if matrix given
  if(is.null(y)) {
    if(!is.matrix(x)) 
      stop("Single argument requires a global cost matrix");
    
    lm <- x;
  } else {
    ## two ts. given TODO handle multivariate
    ## as.numeric handles ts, but not multivariates
    lm <- proxy::dist(x,y,method=dist.method);
  }


  ## Now we have a function
  wfun<-.canonicalizeWindowFunction(window.type);
  

  ## Now we have a step pattern
  dir<-.canonicalizeStepPattern(step.pattern);


  ## shorthand names
  n <- nrow(lm);
  m <- ncol(lm);

  ## perform the computation
  gcm <- globalCostMatrix(lm, step.matrix=dir,
                          window.function=wfun, ...);


  
  jmin <- m;


  ## remember size
  gcm$N <- n;
  gcm$M <- m;

  ## store call
  gcm$call <- match.call();

  ## result: distance (add to existing list gcm?)
  gcm$distance <- gcm$costMatrix[n,jmin];

  ## alignment valid?
  if(is.na(gcm$distance))
    stop("No warping paths exists that is allowed by costraints"); 

  
  if(!distance.only) {
    ## perform the backtrack
    mapping <- backtrack(jmin,gcm);

    ## append to existing list gcm, for now
    ## perhaps replace by attr()
    gcm$index1 <- mapping$index1;
    gcm$index2 <- mapping$index2;
  }


  ## delete sizey intermediate steps 
  if(!keep.internals) {
      gcm$costMatrix<-NULL;
      gcm$directionMatrix<-NULL;
  } else {
      gcm$localCostMatrix <- lm;
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




###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino@gmail.com>           #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: dtw.R 155 2008-06-19 15:12:04Z tonig $
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
         partial=FALSE,
         ... ) {

  lm <- NULL;


  
  ## if matrix given
  if(is.null(y)) {
      if(!is.matrix(x)) 
        stop("Single argument requires a global cost matrix");
    
      lm <- x;
  } else {
      ## two timeseries or vectors given
      ## as.matrix coerces ts or mts to matrices
      x <- as.matrix(x);
      y <- as.matrix(y);
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


  ## remember size
  gcm$N <- n;
  gcm$M <- m;

  ## remember  call
  gcm$call <- match.call();



  ## last column, normalized
  norm <- attr(dir,"norm");
  lastcol <- gcm$costMatrix[n,];

  if(is.na(norm)) {
      # NO-OP
  } else if(norm == "N+M") {
      lastcol <- lastcol/(n+(1:m));
  } else if(norm == "N") {
      lastcol <- lastcol/n;
  } else if(norm == "M") {
      lastcol <- lastcol/m;
  }

  
  ## for complete alignment
  gcm$jmin <- m;

  ## for partial alignment: normalize
  if (partial) {
    if(is.na(norm)) {
        warning("Unknown normalization with partial=TRUE: using N+M");
        lastcol <- lastcol/(n+(1:m));
    }
    gcm$jmin <- which.min(lastcol);
  }

  ## result: distance 
  gcm$distance <- gcm$costMatrix[n,gcm$jmin];

  ## alignment valid?
  if(is.na(gcm$distance)) {
    stop("No warping paths exists that is allowed by costraints"); 
  }
  
  
  ## normalized distance
  if(! is.na(norm)) {
      gcm$normalizedDistance <- lastcol[gcm$jmin];
  } else {
      gcm$normalizedDistance <- NA;
  }

  
  if(!distance.only) {
    ## perform the backtrack
    mapping <- backtrack(gcm);

    gcm$index1 <- mapping$index1;
    gcm$index2 <- mapping$index2;
  }


  ## don't keep internals: delete sizey intermediate steps 
  if(!keep.internals) {
      gcm$costMatrix<-NULL;
      gcm$directionMatrix<-NULL;
  } else {
  ## keep internals: add data
      gcm$localCostMatrix <- lm;
      if(! is.null(y)) {
          gcm$query <- x;
          gcm$template <- y;
      }
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




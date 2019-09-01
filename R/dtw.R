
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
## Frontend stuff, including coercing shorthand types
##




#' Dynamic Time Warp
#' 
#' Compute Dynamic Time Warp and find optimal alignment between two time
#' series.
#' 
#' 
#' The function performs Dynamic Time Warp (DTW) and computes the optimal
#' alignment between two time series `x` and `y`, given as numeric
#' vectors.  The "optimal" alignment minimizes the sum of distances between
#' aligned elements. Lengths of `x` and `y` may differ.
#' 
#' The local distance between elements of `x` (query) and `y`
#' (reference) can be computed in one of the following ways:
#' 
#'  1. if `dist.method` is a string, `x` and `y` are passed to the [proxy::dist()] function in package \pkg{proxy} with the method given; 
#'  2. if `dist.method` is a function of two arguments, it invoked repeatedly on all pairs `x[i],y[j]` to build the local cost matrix; 
#'  3. multivariate time series and arbitrary distance metrics can be handled by supplying a local-distance matrix. Element `[i,j]` of the local-distance matrix is understood as the distance between element `x[i]` and `y[j]`. The distance matrix has therefore `n=length(x)` rows and `m=length(y)` columns (see note below).  
#' 
#' Several common variants of the DTW recursion are supported via the
#' `step.pattern` argument, which defaults to `symmetric2`. Step
#' patterns are commonly used to *locally* constrain the slope of the
#' alignment function. See [stepPattern()] for details.
#' 
#' Windowing enforces a *global* constraint on the envelope of the warping
#' path. It is selected by passing a string or function to the
#' `window.type` argument. Commonly used windows are (abbreviations
#' allowed):
#' 
#'  * `"none"` No windowing (default) 
#'  * `"sakoechiba"` A band around main diagonal 
#'  * `"slantedband"` A band around slanted diagonal 
#'  * `"itakura"` So-called Itakura parallelogram 
#' 
#' `window.type` can also be an user-defined windowing function.  See
#' [dtwWindowingFunctions()] for all available windowing functions,
#' details on user-defined windowing, and a discussion of the (mis)naming of
#' the "Itakura" parallelogram as a global constraint.  Some windowing
#' functions may require parameters, such as the `window.size` argument.
#' 
#' Open-ended alignment, i.e. semi-unconstrained alignment, can be selected via
#' the `open.end` switch.  Open-end DTW computes the alignment which best
#' matches all of the query with a *leading part* of the reference. This
#' is proposed e.g. by Mori (2006), Sakoe (1979) and others. Similarly,
#' open-begin is enabled via `open.begin`; it makes sense when
#' `open.end` is also enabled (subsequence finding). Subsequence
#' alignments are similar e.g. to UE2-1 algorithm by Rabiner (1978) and others.
#' Please find a review in Tormene et al. (2009).
#' 
#' If the warping function is not required, computation can be sped up enabling
#' the `distance.only=TRUE` switch, which skips the backtracking step. The
#' output object will then lack the `index{1,2,1s,2s}` and
#' `stepsTaken` fields.
#' 
#' `is.dtw` tests whether the argument is of class `dtw`.
#' 
#' @aliases is.dtw print.dtw
#' @param x query vector *or* local cost matrix
#' @param y reference vector, unused if `x` given as cost matrix
#' @param dist.method pointwise (local) distance function to use. See
#' [proxy::dist()] in package \pkg{proxy}
#' @param step.pattern a stepPattern object describing the local warping steps
#' allowed with their cost (see [stepPattern()])
#' @param window.type windowing function. Character: "none", "itakura",
#' "sakoechiba", "slantedband", or a function (see details).
#' @param open.begin,open.end perform open-ended alignments
#' @param keep.internals preserve the cumulative cost matrix, inputs, and other
#' internal structures
#' @param distance.only only compute distance (no backtrack, faster)
#' @param d an arbitrary R object
#' @param ... additional arguments, passed to `window.type`
#' @return An object of class `dtw` with 
#' the following items:
#' 
#'  * `distance` the minimum global distance computed, *not* normalized.
#'  * `normalizedDistance` distance computed, *normalized* for path length, if normalization is known for chosen step pattern.
#'  * `N,M` query and reference length
#'  * `call` the function call that created the object
#'  * `index1` matched elements: indices in `x`
#'  * `index2` corresponding mapped indices in `y`
#'  * `stepPattern` the `stepPattern` object used for the computation
#'  * `jmin` last element of reference matched, if `open.end=TRUE`
#'  * `directionMatrix` if `keep.internals=TRUE`, the directions of steps that would be taken at each alignment pair (integers indexing  production rules in the chosen step pattern)
#'  * `stepsTaken` the list of steps taken from the beginning to the end of the alignment (integers indexing chosen step pattern)
#'  * `index1s, index2s` same as `index1/2`, excluding intermediate steps for multi-step patterns like [asymmetricP05()] 
#'  * `costMatrix` if `keep.internals=TRUE`, the cumulative cost matrix
#'  * `query, reference` if `keep.internals=TRUE` and passed as the `x` and `y` arguments, the query and reference timeseries.
#' @note Cost matrices (both input and output) have query elements arranged
#' row-wise (first index), and reference elements column-wise (second index).
#' They print according to the usual convention, with indexes increasing down-
#' and rightwards.  Many DTW papers and tutorials show matrices according to
#' plot-like conventions, i.e.  reference index growing upwards. This may be
#' confusing.
#' @author Toni Giorgino
#' @seealso [dtwDist()], for iterating dtw over a set of timeseries;
#' [dtwWindowingFunctions()], for windowing and global constraints;
#' [stepPattern()], step patterns and local constraints;
#' [plot.dtw()], plot methods for DTW objects.  To generate a local
#' distance matrix, the functions [proxy::dist()] in package
#' \pkg{proxy}, [analogue::distance()] in package \pkg{analogue},
#' [outer()] may come handy.
#' @references
#' 1. Toni Giorgino. *Computing and Visualizing Dynamic Time
#' Warping Alignments in R: The dtw Package.* Journal of Statistical Software,
#' 31(7), 1-24. <http://www.jstatsoft.org/v31/i07/>
#' 2. Tormene, P.;
#' Giorgino, T.; Quaglini, S. & Stefanelli, M. *Matching incomplete time
#' series with dynamic time warping: an algorithm and an application to
#' post-stroke rehabilitation.* Artif Intell Med, 2009, 45, 11-34.
#' <http://dx.doi.org/10.1016/j.artmed.2008.11.007>
#' 3. Sakoe, H.;
#' Chiba, S., *Dynamic programming algorithm optimization for spoken word
#' recognition,* Acoustics, Speech, and Signal Processing,
#' IEEE Transactions on , vol.26, no.1, pp. 43-49, Feb 1978.
#' <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055> 
#' 4. Mori, A.; Uchida, S.; Kurazume, R.; Taniguchi, R.; Hasegawa, T. & Sakoe, H.
#' *Early Recognition and Prediction of Gestures* Proc. 18th International
#' Conference on Pattern Recognition ICPR 2006, 2006, 3, 560-563
#' 5. Sakoe,
#' H. *Two-level DP-matching--A dynamic programming-based pattern matching
#' algorithm for connected word recognition* Acoustics, Speech, and Signal
#' Processing, IEEE
#' Transactions on, 1979, 27, 588-595
#' 6. Rabiner L, Rosenberg A, Levinson
#' S (1978). *Considerations in dynamic time warping algorithms for
#' discrete word recognition.* IEEE Trans. Acoust., Speech, Signal Process.,
#' 26(6), 575-582. ISSN 0096-3518.
#' 7. Muller M. *Dynamic Time
#' Warping* in *Information Retrieval for Music and Motion*. Springer
#' Berlin Heidelberg; 2007. p. 69-84.
#' <http://link.springer.com/chapter/10.1007/978-3-540-74048-3_4>
#' @keywords ts
#' @concept Dynamic Time Warp
#' @concept Dynamic programming
#' @concept Align timeseries
#' @concept Minimum cumulative cost
#' @concept Distance
#' @examples
#' 
#' 
#' ## A noisy sine wave as query
#' idx<-seq(0,6.28,len=100);
#' query<-sin(idx)+runif(100)/10;
#' 
#' ## A cosine is for reference; sin and cos are offset by 25 samples
#' reference<-cos(idx)
#' plot(reference); lines(query,col="blue");
#' 
#' ## Find the best match
#' alignment<-dtw(query,reference);
#' 
#' 
#' ## Display the mapping, AKA warping function - may be multiple-valued
#' ## Equivalent to: plot(alignment,type="alignment")
#' plot(alignment$index1,alignment$index2,main="Warping function");
#' 
#' ## Confirm: 25 samples off-diagonal alignment
#' lines(1:100-25,col="red")
#' 
#' 
#' 
#' 
#' #########
#' ##
#' ## Partial alignments are allowed.
#' ##
#' 
#' alignmentOBE <-
#'   dtw(query[44:88],reference,
#'       keep=TRUE,step=asymmetric,
#'       open.end=TRUE,open.begin=TRUE);
#' plot(alignmentOBE,type="two",off=1);
#' 
#' 
#' #########
#' ##
#' ## Subsetting allows warping and unwarping of
#' ## timeseries according to the warping curve. 
#' ## See first example below.
#' ##
#' 
#' ## Most useful: plot the warped query along with reference 
#' plot(reference)
#' lines(query[alignment$index1]~alignment$index2,col="blue")
#' 
#' ## Plot the (unwarped) query and the inverse-warped reference
#' plot(query,type="l",col="blue")
#' points(reference[alignment$index2]~alignment$index1)
#' 
#' 
#' 
#' #########
#' ##
#' ## Contour plots of the cumulative cost matrix
#' ##    similar to: plot(alignment,type="density") or
#' ##                dtwPlotDensity(alignment)
#' ## See more plots in ?plot.dtw 
#' ##
#' 
#' ## keep = TRUE so we can look into the cost matrix
#' 
#' alignment<-dtw(query,reference,keep=TRUE);
#' 
#' contour(alignment$costMatrix,col=terrain.colors(100),x=1:100,y=1:100,
#' 	xlab="Query (noisy sine)",ylab="Reference (cosine)");
#' 
#' lines(alignment$index1,alignment$index2,col="red",lwd=2);
#' 
#' 
#' 
#' 
#' #########
#' ##
#' ## An hand-checkable example
#' ##
#' 
#' ldist<-matrix(1,nrow=6,ncol=6);  # Matrix of ones
#' ldist[2,]<-0; ldist[,5]<-0;      # Mark a clear path of zeroes
#' ldist[2,5]<-.01;		 # Forcely cut the corner
#' 
#' ds<-dtw(ldist);			 # DTW with user-supplied local
#'                                  #   cost matrix
#' da<-dtw(ldist,step=asymmetric);	 # Also compute the asymmetric 
#' plot(ds$index1,ds$index2,pch=3); # Symmetric: alignment follows
#'                                  #   the low-distance marked path
#' points(da$index1,da$index2,col="red");  # Asymmetric: visiting
#'                                         #   1 is required twice
#' 
#' ds$distance;
#' da$distance;
#' 
#' 
#' 
#' 
#' 
#' @export dtw
`dtw` <-
function(x, y=NULL,
         dist.method="Euclidean",
         step.pattern=symmetric2,
         window.type="none",
         keep.internals=FALSE,
         distance.only=FALSE,
         open.end=FALSE,
         open.begin=FALSE,
         ... ) {

  ## The local cost matrix  
  lm <- NULL;

  ## if matrix given
  if(is.null(y)) {
      if(!missing(dist.method))
          stop("Argument dist.method does not make sense with a local cost matrix")
      if(!is.matrix(x)) 
          stop("Single argument requires a pre-computed local cost matrix");
      lm <- x;
  } else {
      ## two timeseries or vectors given
      ## as.matrix coerces ts or mts to matrices
      x <- as.matrix(x);
      y <- as.matrix(y);
      if( (ncol(x)==1 || ncol(y)==1) && !missing(dist.method) )
          warning("Argument dist.method is only useful with multivariate timeseries")
      if(!is.character(dist.method)) 
          stop("dist.method should be a method name supported by proxy::dist()");
      lm <- proxy::dist(x,y,method=dist.method);
  } 
  

  ## Now we have a function
  wfun<-.canonicalizeWindowFunction(window.type);
  

  ## Now we have a step pattern
  dir<-step.pattern;
  norm <- attr(dir,"norm");


  ## Warn for obsolete constructs
  if(! is.null(list(...)$partial) ) {
    warning("Argument `partial' is obsolete. Use `open.end' instead");
    open.end <- TRUE;
  }



  ## shorthand names
  n <- nrow(lm);
  m <- ncol(lm);

  
  ## For open-begin alignment:
  if (open.begin) {
    
    ##  ensure proper normalization
    if(is.na(norm) || norm != "N") {
      stop("Open-begin requires step patterns with 'N' normalization (e.g. asymmetric, or R-J types (c)). See papers in citation().");
    }

    ## prepend a null row
    lm <- rbind(0,lm);
    np <- n+1;

    ##  pre-initialize elements in the cumulative cost matrix
    precm <- matrix(NA,nrow=np,ncol=m);
    precm[1,] <- 0;

  } else {
    precm <- NULL;
    np <- n;
  }

  
  ## perform the computation
  gcm <- globalCostMatrix(lm, step.matrix=dir,
                          window.function=wfun,
                          seed=precm, ...);


  ## remember size
  gcm$N <- n;
  gcm$M <- m;

  ## remember  call
  gcm$call <- match.call();
  gcm$openEnd <- open.end;
  gcm$openBegin <- open.begin;
  gcm$windowFunction <- wfun;

  ## last row (misnamed), normalized
  lastcol <- gcm$costMatrix[np,];

  if(is.na(norm)) {
      # NO-OP
  } else if(norm == "N+M") {
      lastcol <- lastcol/(n+(1:m));
  } else if(norm == "N") {
      lastcol <- lastcol/n;
  } else if(norm == "M") {
      lastcol <- lastcol/(1:m);
  }

  
  ## for complete alignment
  gcm$jmin <- m;

  
  ## for open-end alignment: normalize
  if (open.end) {
    if(is.na(norm)) {
      stop("Open-end alignments require normalizable step patterns");
    }
    gcm$jmin <- which.min(lastcol);
  }

  ## result: distance 
  gcm$distance <- gcm$costMatrix[np,gcm$jmin];

  ## alignment valid?
  if(is.na(gcm$distance)) {
    stop("No warping path exists that is allowed by costraints"); 
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
    gcm <- c(gcm,mapping);    ## Add the properties to gcm
  }


  ## open-begin: discard first elements
  if(open.begin) {
    gcm$index1 <- gcm$index1[-1]-1;
    gcm$index1s <- gcm$index1s[-1]-1;
    gcm$index2 <- gcm$index2[-1];
    gcm$index2s <- gcm$index2s[-1];
    lm <- lm[-1,];
    gcm$costMatrix <- gcm$costMatrix[-1,];
    gcm$directionMatrix <- gcm$directionMatrix[-1,];
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
          gcm$reference <- y;
      }
  }


  ## if a dtw object is to be sponsored:
  class(gcm) <- "dtw";
  return(gcm);
}


##############################
## OO class check
#' @rdname  dtw  
#' @export
is.dtw <- function(d) {
    return(inherits(d,"dtw"));
}



##############################
## OO print method
#' @rdname dtw 
#' @export
print.dtw <- function(x,...) {
  head <- "DTW alignment object\n";
  size <- sprintf("Alignment size (query x reference): %d x %d\n",x$N,x$M);
  call <- sprintf("Call: %s\n",deparse(x$call));
  cat(head,size,call);
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




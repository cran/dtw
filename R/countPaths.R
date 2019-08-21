###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino,gmail.com>               #
#       Istituto di Neuroscienze (IN-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id$
#                                                             #
###############################################################






#' Count the number of possible warping paths
#' 
#' Count the number of warping paths compatible with the constraints.
#' 
#' Count how many possible warping paths exist in the alignment problem passed
#' as an argument. The object passed as an argument is used to look up the
#' problem parameters such as the used step pattern, windowing, open ends, and
#' so on. The actual alignment is ignored.
#' 
#' Note that the number of paths grows exponentially with problems size. The
#' result may be approximate when windowing functions are used.
#' 
#' If \code{debug} is \code{TRUE}, a matrix used for the computation is
#' returned instead of the final result.
#' 
#' @param d an object of class \code{dtw}
#' @param debug return an intermediate result
#' @return The number of paths.
#' @author Toni Giorgino
#' @keywords ts
#' @examples
#' 
#'   ds<-dtw(1:7+2,1:8,keep=TRUE,step=asymmetric);
#'   countPaths(ds)
#'   ## Result: 126
#' 
#' @export countPaths
countPaths <-
function(d,debug=FALSE) {
  
  N<-d$N;
  M<-d$M;
  m<-matrix(NA,nrow=N,ncol=M);
  
  if(d$openBegin) {
    m[1,]<-  1;
  } else {
    m[1,1]<-1;
  }
  
  # Some help functions
  dir<-d$stepPattern;
  npats <- attr(dir,"npat");
  nsteps <- dim(dir)[1];          # number of individual steps (counting all patterns)
  deltas <- .mkDirDeltas(dir);     # Cache the total step for each pattern

  wf <- d$windowFunction;
  
  for (ii in 1:N) {
    for (jj in 1:M) {

      # do nothing if already computed
      if( is.finite(m[ii,jj]) ) { next }

      # set to zero if outside window. Ugly but necessary to recover
      # optional arguments such as window.size
      if( do.call(d$windowFunction,
                  c(list(iw=ii,jw=jj),
                  as.list(d$call)))==FALSE ) {
        m[ii,jj]<-0; 
        next;
      }

      np<-0;
      for (k in 1:npats) {
        ni<-ii-deltas[k,1];
        nj<-jj-deltas[k,2];

        if(ni>=1 && nj>=1 ) {
          np <- np+m[ni,nj];
        }
      }
      m[ii,jj] <- np
    }
  }

  if(debug) {
    return(m)
  } 
  
  if(d$openEnd) {
    return(sum(m[N,]));
  } else {
    return(m[N,M]);
  }
  
}

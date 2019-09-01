
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






#' Count the number of warping paths consistent with the constraints.
#' 
#' Count how many possible warping paths exist in the alignment problem passed
#' as an argument. The object passed as an argument is used to look up the
#' problem parameters such as the used step pattern, windowing, open ends, and
#' so on. The actual alignment is ignored.
#' 
#' Note that the number of paths grows exponentially with problems size. The
#' result may be approximate when windowing functions are used.
#' 
#' If `debug=TRUE`, a matrix used for the computation is
#' returned instead of the final result.
#' 
#' @param d an object of class `dtw`
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


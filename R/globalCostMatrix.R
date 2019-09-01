
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


########################################
## Compute the cost matrix from a local distance matrix

## Wrapper to the native function

#' @useDynLib dtw, .registration = TRUE, .fixes="C_" 



`globalCostMatrix` <-
function(lm,
         step.matrix=symmetric1,
         window.function=noWindow,
         native=TRUE,
         seed=NULL,
         ...) {


  ## sanity check - be extra cautions w/ binary
  if (!is.stepPattern(step.matrix))
    stop("step.matrix is no stepMatrix object");



  # i = 1 .. n in query sequence, on first index, ie rows
  # j = 1 .. m on reference sequence, on second index, ie columns
  #   Note:  reference is usually drawn vertically, up-wise

  n <- nrow(lm);
  m <- ncol(lm);


  # number of individual steps (counting all patterns)
  nsteps<-dim(step.matrix)[1];


  # clear the cost and step matrix
  # these will be the outputs of the binary
  # for  cm use  seed if given
  if(!is.null(seed)) {
    cm <- seed;
  } else {
    cm <- matrix(NA,nrow=n,ncol=m);
    cm[1,1] <- lm[1,1];                 # Questionable.
  }

  sm <- matrix(NA,nrow=n,ncol=m);

  #if(is.loaded("computeCM_Call") && native){ # -- not working any more?
  if(native){
        ## precompute windowing
        wm <- matrix(FALSE,nrow=n,ncol=m);
        wm[window.function(row(wm),col(wm),
                       query.size=n, reference.size=m,
                       ...)]<-TRUE;
        
        storage.mode(wm) <- "logical";
        storage.mode(lm) <- "double";
        storage.mode(cm) <- "double";
        storage.mode(step.matrix) <- "double";
        out <- .Call(C_computeCM_Call,
                   wm,lm,cm,step.matrix);
        
  } else {
    ####################
    ## INTERPRETED PURE-R IMPLEMENTATION
    warning("Native dtw implementation not available: using (slow) interpreted fallback");
                                        # now walk through the matrix, column-wise and row-wise,
                                        # and recursively compute the accumulated distance. Unreachable
                                        # elements are handled via NAs (removed)
    dir <- step.matrix;
    npats <- attr(dir,"npat");
    for (j in 1:m) {
      for (i in 1:n) {
        ## It is ok to window on the arrival point (?)
        if(!window.function(i,j, query.size=n, reference.size=m, ...)) { next; }

        ## Skip if already initialized
        if(!is.na(cm[i,j])) { next; }

        clist<-numeric(npats)+NA;
        for (s in 1:nsteps) {
          ## current pattern
          p<-dir[s,1];
          ## ii,jj is the cell from where potentially we could
          ## have come from. 
          ii<-i-dir[s,2];                 # previous step in inp
          jj<-j-dir[s,3];                 # previous step in tpl
          if(ii>=1 && jj>=1) {            # element exists?
            cc<-dir[s,4];                 # step penalty
            if(cc == -1) {                #  -1? cumulative cost:
              clist[p]<-cm[ii,jj];	#  there must be exactly 1 per pattern
            } else {			#  a cost for 
              clist[p]<-clist[p]+cc*lm[ii,jj];
            }
          }
        }


        ## no NAs in clist at this point BUT clist can be empty
        ## store in cost matrix
        minc<-which.min(clist);           # pick the least cost
        if(length(minc)>0) {          	# false if clist has all NAs
          cm[i,j]<-clist[minc];
          sm[i,j]<-minc;			# remember the pattern picked
        }
      }
    }
    out <- list(costMatrix=cm,directionMatrix=sm);
  }

  ## END PURE-R IMPLEMENTATION
  ####################

  ## At this point out$cmo and out$smo should be set
  out$stepPattern <- step.matrix;
  return(out);
}



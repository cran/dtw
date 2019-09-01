
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
## Backtrack the steps taken - internal

`backtrack` <-
function(gcm) {

  dir<-gcm$stepPattern;
  npat <- attr(dir,"npat");

  n <- nrow(gcm$costMatrix);
  m <- ncol(gcm$costMatrix);

  i <- n;
  j <- gcm$jmin;


  ## drop rows with (0,0) deltas 
  nullrows <- dir[,2]==0 & dir[,3] ==0 ;
  tmp <- dir[!nullrows,,drop=FALSE];

  ## Pre-compute steps
  stepsCache <- list();  
  for(k in 1:npat) {
    stepsCache[[k]] <- .extractpattern(tmp,k);
  }


  ## mapping lists
  iis<-ii<-c(i);
  jjs<-jj<-c(j);
  ss<-NULL;

  repeat {
    ## cross fingers for termination
    if(i==1 && j==1) {
      break;
    }

    ## direction taken
    s<-gcm$directionMatrix[i,j];

    if(is.na(s)) {
      break;
    }

    ## undo the steps
    ss<-c(s,ss);
    steps<-stepsCache[[s]];
    ns<-nrow(steps);

    ## In some rare cases (eg symmetricP0), ns will be 1
    ## R indexing rules make k==0 a no-op anyway
    for(k in 1:ns) {
	## take note of current cell, prepending to mapping lists
	ii <- c(i-steps[k,1],ii);
	jj <- c(j-steps[k,2],jj);
	## All sub-steps are visited & appended; we have dropped (0,0) deltas
    }

    ## And don't forget where we arrived to
    i<-i-steps[ns,1];
    j<-j-steps[ns,2];

    ## iis & jjs are backtracking without the intermediate steps
    iis<-c(i,iis)
    jjs<-c(j,jjs)
  }


  out<-list(index1=ii,index2=jj,
            index1s=iis, index2s=jjs,
            stepsTaken=ss);

  return(out);
}


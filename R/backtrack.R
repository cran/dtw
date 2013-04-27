###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino,gmail.com>               #
#       Istituto di Ingegneria Biomedica (ISIB-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id: backtrack.R 286 2013-04-27 19:10:00Z tonig $
#                                                             #
###############################################################


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
  ii<-c(i);
  jj<-c(j);
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
  }


  out<-list(index1=ii,index2=jj,stepsTaken=ss);

  return(out);
}

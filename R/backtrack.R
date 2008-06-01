###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino@gmail.com>               #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: backtrack.R 126 2008-05-29 12:44:37Z tonig $
#                                                             #
###############################################################


########################################
## Backtrack the steps taken - internal

`backtrack` <-
function(gcm) {

  dir<-gcm$stepPatterns;
  npat <- attr(dir,"npat");

  n <- nrow(gcm$costMatrix);
  m <- ncol(gcm$costMatrix);

  i <- n;
  j <- gcm$jmin;



  ## drop rows with (0,0) deltas 
  nullrows <- dir[,2]==0 & dir[,3] ==0 ;
  tmp <- dir[!nullrows,];

  ## Pre-compute steps
  stepsCache <- list();  
  for(k in 1:npat) {
    stepsCache[[k]] <- .extractpattern(tmp,k);
  }

  
  ## mapping lists
  ii<-c(i);
  jj<-c(j);



  repeat {
    ## cross fingers for termination
    if(i==1 && j==1) {
      break;
    }

    ## direction taken
    s<-gcm$directionMatrix[i,j];

    ## undo the steps

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


  out<-list();
  out$index1<-ii;
  out$index2<-jj;

  return(out);
}

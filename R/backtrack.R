###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino@gmail.com>               #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: backtrack.R 31 2007-12-07 10:46:09Z tonig $
#                                                             #
###############################################################


########################################
## Backtrack the steps taken - internal

`backtrack` <-
function(jmin, gcm) {

  ## mapping lists
  ii<-c();
  jj<-c();

  dir<-gcm$stepPatterns;

  n <- dim(gcm$costMatrix)[1];
  m <- dim(gcm$costMatrix)[2];

  i <- n;
  j <- jmin;

  repeat {
    ## cross fingers for termination
    if(i==1 && j==1) {
      break;
    }

    ## direction taken
    s<-gcm$directionMatrix[i,j];

    ## undo the steps

    steps<-.extractpattern(dir,s);
    ns<-dim(steps)[1];

    for(k in 1:(ns-1)) {
	## take note of current cell, prepending to mapping lists
	ii <- c(i-steps[k,1],ii);
	jj <- c(j-steps[k,2],jj);
	## The final cell is not appended; it will be done 
	## processing step 0,0 at the next iteration
    }

    ## And don't forget where we arrived
    i<-i-steps[ns,1];
    j<-j-steps[ns,2];
  }

  ## The very final point of the backtrack is not otherwise appended 
  ii<-c(1,ii);
  jj<-c(1,jj);

  out<-list();
  out$index1<-ii;
  out$index2<-jj;

  return(out);
}




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
    ## take note of current cell, prepending to mapping lists
    ii <- c(i,ii);
    jj <- c(j,jj);

    ## cross fingers for termination
    if(i==1 && j==1) {
      break;
    }

    ## direction taken
    s<-gcm$directionMatrix[i,j];

    ## undo the step
    i<-i-dir[s,1];
    j<-j-dir[s,2];

  }

  out<-list();
  out$index1<-ii;
  out$index2<-jj;

  return(out);
}


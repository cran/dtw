
########################################
## Compute the cost matrix from a local distance matrix
## $Id:$





`globalCostMatrix` <-
function(lm,
         step.pattern="s",
         window.type="none",
         window.size=10) {


  ## TODO. Coerce lm to matrix?
  
  
  wfun<-window.type;

  if(!is.function(wfun)) {
    wt<-pmatch(window.type,c("sakoechiba","itakura","none"));    
    if(is.na(wt)) {
      wfun<-noWindow;
      warning("Unsupported argument for window.type: using no windowing");
    } else if(wt==1) {                  #SAKOE band
      wfun<-sakoeChibaWindow;
    } else if(wt==2) {                  #ITAKURA pg.
      wfun<-itakuraWindow;
    } else if(wt==3) {                  #no windowing
      wfun<-noWindow;
    } else {
      wfun<-noWindow;
      warning("Should not happen on window.type: using no windowing");
    }
  }
  

  # Step patterns: \delta i, \delta j, cost
  #                                      ...
  # all deltas MUST be positive (otherwise we violate monotonicity)
  #       #TODO pmatch
  
  dir<-step.pattern;        # custom step pattern
  if(dir == "s") { # symmetric DTW
    dir<-matrix(c(0,1,1,
                  1,0,1,
                  1,1,1), nrow=3, ncol=3,byrow=TRUE);
  } else if(dir == "a") { # asymmetric DTW
    dir<-matrix(c(1,0,1,
                  1,1,1,
                  1,2,1), nrow=3, ncol=3,byrow=TRUE);
  }

  # i = 1 .. n in query sequence, on first index, ie rows
  # j = 1 .. m on template sequence, on second index, ie columns
  #   Note:  template is usually drawn vertically, up-wise

  query.size    <- n <- dim(lm)[1];
  template.size <- m <- dim(lm)[2];


  # number of step patterns allowed
  nsteps<-dim(dir)[1];

  # clear the cost and step matrix
  cm <- matrix(NA,nrow=n,ncol=m);
  sm <- matrix(NA,nrow=n,ncol=m);

  # initializer
  cm[1,1] <- lm[1,1];


  # now walk through the matrix, column-wise and row-wise,
  # and recursively compute the accumulated distance. Unreachable
  # elements are handled via NAs (removed)

  for (j in 1:m) {
    for (i in 1:n) {
      d<-lm[i,j];
      clist<-c();
      slist<-c();
      for (s in 1:nsteps) {
        ## ii,jj is the cell from where potentially we could
        ## come from. TODO add windowing here.
        ii<-i-dir[s,1];                 # previous step in inp
        jj<-j-dir[s,2];                 # previous step in tpl
        if(ii>=1 && jj>=1) {            # element exists?
          cc<-  dir[s,3];                 # step penalty
          clist<-c(clist,cm[ii,jj]+cc*d); #build a cost list
          slist<-c(slist,s);              #remember whence we came
        }
      }


      ## there could no NAs at this point

      # store in cost matrix
      if(length(which.min(clist))>0) {          # odd way to exclude NANs
        minc<-which.min(clist);                 # pick the least cost
        cm[i,j]<-clist[minc];
        sm[i,j]<-slist[minc];
      }
    }
  }


  out<-list();
  out$costMatrix<-cm;                   # to get distance
  out$directionMatrix<-sm;              # to backtrack
  out$stepPatterns<-dir;                # to backtrack

  return(out);
}


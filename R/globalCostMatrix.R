###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino@gmail.com>               #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: globalCostMatrix.R 51 2007-12-11 10:59:08Z tonig $
#                                                             #
###############################################################


########################################
## Compute the cost matrix from a local distance matrix



## We assume that all arguments are expanded from
## char shortands. This includes:
##  step.matrix - should be a stepMatrix object
##  window.function - should be a function





`globalCostMatrix` <-
function(lm,
         step.matrix=symmetric1,
         window.function=noWindow,
         ...) {


  # a shorthand
  dir <- step.matrix;

  # i = 1 .. n in query sequence, on first index, ie rows
  # j = 1 .. m on template sequence, on second index, ie columns
  #   Note:  template is usually drawn vertically, up-wise

  query.size    <- n <- dim(lm)[1];
  template.size <- m <- dim(lm)[2];


  # number of individual steps (counting all patterns)
  nsteps<-dim(dir)[1];

  # number of step patterns defined
  npats<-max(dir[,1]); 

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
      ## It is ok to window on the arrival point (?)
      if(!window.function(i,j, query.size, template.size, ...)) { next; }

      clist<-numeric(npats)+NA;
      for (s in 1:nsteps) {
	## current pattern
	p<-dir[s,1];
        ## ii,jj is the cell from where potentially we could
        ## have come from. 
        ii<-i-dir[s,2];                 # previous step in inp
        jj<-j-dir[s,3];                 # previous step in tpl
        if(ii>=1 && jj>=1) {            # element exists?
          cc<-  dir[s,4];               # step penalty
	  if(cc == -1) {		#  -1? cumulative cost:
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


  out<-list();
  out$costMatrix<-cm;                   # to get distance
  out$directionMatrix<-sm;              # to backtrack
  out$stepPatterns<-dir;                # to backtrack

  return(out);
}


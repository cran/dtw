###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino@gmail.com>               #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: globalCostNative.R 51 2007-12-11 10:59:08Z tonig $
#                                                             #
###############################################################


########################################
## Compute the cost matrix from a local distance matrix

## Wrapper to the native function





`globalCostNative` <-
function(lm,
         step.matrix=symmetric1,
         window.function=noWindow,
         ...) {


  ## sanity check - be extra cautions w/ binary
  if (!is.stepPattern(step.matrix))
    stop("step.matrix is no stepMatrix object");



  # i = 1 .. n in query sequence, on first index, ie rows
  # j = 1 .. m on template sequence, on second index, ie columns
  #   Note:  template is usually drawn vertically, up-wise

  n <- dim(lm)[1];
  m <- dim(lm)[2];


  # number of individual steps (counting all patterns)
  nsteps<-dim(step.matrix)[1];


  # clear the cost and step matrix
  # these will be the outputs of the binary
  cm <- matrix(NA,nrow=n,ncol=m);
  sm <- matrix(NA,nrow=n,ncol=m);


  ## precompute windowing
  wm <- matrix(FALSE,nrow=n,ncol=m);
  wm[window.function(row(wm),col(wm),
	query.size=n, template.size=m,
	...)]<-TRUE;
  
  
  # initializer
  cm[1,1] <- lm[1,1];

  
  tmp<-.C("computeCM",NAOK=TRUE,PACKAGE="dtw",
     as.integer(dim(cm)),               # s
     as.logical(wm),                    #
     as.double(lm),
     as.integer(nsteps),
     as.double(step.matrix),
     cmo=as.double(cm),                     # OUT
     smo=as.integer(sm));                   # OUT

  cm<-matrix(tmp$cmo,nrow=n,ncol=m);
  sm<-matrix(tmp$smo,nrow=n,ncol=m);


  out<-list();
  out$costMatrix<-cm;                   # to get distance
  out$directionMatrix<-sm;              # to backtrack
  out$stepPatterns<-step.matrix;        # to backtrack

  return(out);
}


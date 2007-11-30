###############################################################
#                                                             #
#   Author: Toni Giorgino <toni.giorgino@gmail.com>           #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: dtw.R 10 2007-12-03 19:17:59Z tonig $
#                                                             #
###############################################################


## For pre-defined step patterns see below.


#############################
## Methods for accessing and creating step.patterns

stepPattern <- function(v) {
  if(!is.vector(v)) {
    stop("stepPattern creation only supported from vectors");
  }

  obj<-matrix(v,ncol=4,byrow=TRUE);
  class(obj)<-"stepPattern";
  return(obj);
}


is.stepPattern <- function(x) {
  return(inherits(x,"stepPattern"));
}

  




## pretty-print the matrix meaning,
## so it will not be as write-only as now

print.stepPattern <-function(x,...) {

  step.pattern<-x;                      # for clarity
  np<-max(step.pattern[,1]);            #no. of patterns

  head<-"g[i,j] = min(\n";
  body<-"";

  ## cycle over available step patterns 
  for(p in 1:np) {
    steps<-.extractpattern(step.pattern,p);
    ns<-dim(steps)[1];

    ## restore row order
    steps<-steps[ns:1,];

    ## cycle over steps s in the current pattern p
    for(s in 1:ns) {
      di<-steps[s,1];                   # delta in query
      dj<-steps[s,2];                   # delta in templ
      cc<-steps[s,3];                   # step cost multiplier

      ## make pretty-printable negative increments
      dis<-ifelse(di==0,"",-di);        # 4 -> -4; 0 -> .
      djs<-ifelse(dj==0,"",-dj);        #  0 maps to empty string

      ## cell origin, as coordinate pair
      dijs<-sprintf("i%2s,j%2s",dis,djs);

      if(cc==-1) {                      # g
        gs<-sprintf("   g[%s]",dijs);
        body<-paste(body,gs);
      } else {
        ## prettyprint step cost multiplier in ccs:  1 -> .; 2 -> 2 *
        ccs<-ifelse(cc==1,"",sprintf(" %d *",cc)); 
        ds<-sprintf("+%s d[%s]",ccs,dijs);
        body<-paste(body,ds);
      }

    }

    body<-paste(body,",\n",s="");
  }

  tail<-")\n\n";

  rv<-paste(head,body,tail);

  cat("Step pattern recursion:\n");
  cat(rv);
  
}



## TODO: sanity check on the step pattern definition

.checkpattern <- function(sp) {
  ## must have 4 x n elements
  ## all integers
  ## first column in ascending order from 1, no missing steps
  ## 2nd, 3rd row non-negative
  ## 4th: first time for each step is -1
}
  

                          


## Extract rows belonging to pattern no. sn
## with first element stripped
## in reverse order

.extractpattern <- function(sp,sn) {
	sbs<-sp[,1]==sn;	# pick only rows beginning by sn
	spl<-sp[sbs,-1];	# of those: take only column Di, Dj, cost
                                # (drop pattern no. column)
        

        ## make sure it stays a matrix
        spl <- matrix(spl,ncol=3);
        
	nr<-dim(spl)[1];	# how many are left
        spl<-spl[nr:1,];	# invert row order

        ## make sure it stays a matrix
        spl <- matrix(spl,ncol=3);

	return(spl);
}


##################################################
##################################################


##
## Various step patterns, defined as internal variables
##

## Some knowledge of DP is required to modify this file.
## Step patterns taken from Sakoe, cited in documentation

## First column: enumerates step patterns.
## Second   	 step in query index
## Third	 step in template index
## Fourth	 weight if positive, or -1 if starting point

## normalization: no
symmetric1 <- stepPattern(c(
                            1,0,1,-1,
                            1,0,0,1,
                            2,1,1,-1,
                            2,0,0,1,
                            3,1,0,-1,
                            3,0,0,1
                            ));

## normalization: N+M
symmetric2 <- stepPattern(c(
                            1,0,1,-1,
                            1,0,0,1,
                            2,1,1,-1,
                            2,0,0,2,
                            3,1,0,-1,
                            3,0,0,1
                            ));


## this one works
## normalization: N
asymmetric <-  stepPattern(c(
                             1,1,0,-1,
                             1,0,0,1,
                             2,1,1,-1,
                             2,0,0,1,
                             3,1,2,-1,
                             3,0,0,1
                           ));


## normalization: N
asymmetricItakura <-  stepPattern(c(
                        1, 1, 2, -1,
			1, 0, 0, 1,
			2, 1, 1, -1,
			2, 0, 0, 1,
			3, 2, 1, -1,
			3, 1, 0, 1,
			3, 0, 0, 1,
			4, 2, 2, -1,
			4, 1, 0, 1,
			4, 0, 0, 1
                       ));




################################

## according to sakoe page 47
## but I'm not very positive about it

.asymmetricSakoe <- stepPattern(c(
                                  1,0,1,-1,
                                  2,1,1,-1,
                                  2,0,0,1,
                                  3,1,0,-1,
                                  3,0,0,1
                                ));




#############################
## Slope-limited versions
##
## Taken from Table I, page 47 of "Dynamic programming algorithm
## optimization for spoken word recognition," Acoustics, Speech, and
## Signal Processing, vol.26, no.1, pp. 43-49, Feb 1978 URL:
## http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055



## Implementation of Sakoe's P=1, Symmetric algorithm 

symmetricP1 <- stepPattern(c(
                              1,1,2,-1,	# First branch: g(i-1,j-2)+
                              1,0,1,2,	#            + 2d(i  ,j-1)  
                              1,0,0,1,	#            +  d(i  ,j)
                              2,1,1,-1,	# Second branch: g(i-1,j-1)+
                              2,0,0,2,	#              +2d(i,  j)
                              3,2,1,-1,	# Third branch: g(i-2,j-1)+
                              3,1,0,2,	#            + 2d(i-1,j)
                              3,0,0,1	#            +  d(  i,j)		
                        ));





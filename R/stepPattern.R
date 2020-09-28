
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


## For pre-defined step patterns see below.


#############################
## Methods for accessing and creating step.patterns
## TODO: validate norm





#' Step patterns for DTW
#' 
#' A `stepPattern` object lists the transitions allowed while searching
#' for the minimum-distance path.  DTW variants are implemented by passing one
#' of the objects described in this page to the `stepPattern` argument of
#' the [dtw()] call.
#' 
#' A step pattern characterizes the matching model and slope constraint
#' specific of a DTW variant. They also known as local- or slope-constraints,
#' transition types, production or recursion rules (GiorginoJSS).
#' 
#' **Pre-defined step patterns**
#' 
#' ```
#'    ## Well-known step patterns
#'    symmetric1
#'    symmetric2
#'    asymmetric
#'    
#'    ## Step patterns classified according to Rabiner-Juang (Rabiner1993)
#'    rabinerJuangStepPattern(type,slope.weighting="d",smoothed=FALSE)
#'    
#'    ## Slope-constrained step patterns from Sakoe-Chiba (Sakoe1978)
#'    symmetricP0;  asymmetricP0
#'    symmetricP05; asymmetricP05
#'    symmetricP1;  asymmetricP1
#'    symmetricP2;  asymmetricP2
#'    
#'    ## Step patterns classified according to Rabiner-Myers (Myers1980)
#'    typeIa;   typeIb;   typeIc;   typeId;
#'    typeIas;  typeIbs;  typeIcs;  typeIds;  # smoothed
#'    typeIIa;  typeIIb;  typeIIc;  typeIId;
#'    typeIIIc; typeIVc;
#'    
#'    ## Miscellaneous
#'    mori2006;
#'    rigid;
#'  ```
#' 
#'  
#' A variety of classification schemes have been proposed for step patterns, including
#' Sakoe-Chiba (Sakoe1978); Rabiner-Juang (Rabiner1993); and Rabiner-Myers
#' (Myers1980).  The `dtw` package implements all of the transition types
#' found in those papers, with the exception of Itakura's and
#' Velichko-Zagoruyko's steps, which require subtly different algorithms (this
#' may be rectified in the future). Itakura recursion is almost, but not quite,
#' equivalent to `typeIIIc`.
#' 
#' For convenience, we shall review pre-defined step patterns grouped by
#' classification. Note that the same pattern may be listed under different
#' names. Refer to paper (GiorginoJSS) for full details.
#' 
#' **1. Well-known step patterns**
#' 
#' Common DTW implementations are based on one of the following transition
#' types.
#' 
#' `symmetric2` is the normalizable, symmetric, with no local slope
#' constraints.  Since one diagonal step costs as much as the two equivalent
#' steps along the sides, it can be normalized dividing by `N+M`
#' (query+reference lengths). It is widely used and the default.
#' 
#' `asymmetric` is asymmetric, slope constrained between 0 and 2. Matches
#' each element of the query time series exactly once, so the warping path
#' `index2~index1` is guaranteed to be single-valued.  Normalized by
#' `N` (length of query).
#' 
#' `symmetric1` (or White-Neely) is quasi-symmetric, no local constraint,
#' non-normalizable. It is biased in favor of oblique steps.
#' 
#' **2. The Rabiner-Juang set**
#' 
#' A comprehensive table of step patterns is proposed in Rabiner-Juang's book
#' (Rabiner1993), tab. 4.5.  All of them can be constructed through the
#' `rabinerJuangStepPattern(type,slope.weighting,smoothed)` function.
#' 
#' The classification foresees seven families, labelled with Roman numerals
#' I-VII; here, they are selected through the integer argument `type`.
#' Each family has four slope weighting sub-types, named in sec. 4.7.2.5 as
#' "Type (a)" to "Type (d)"; they are selected passing a character argument
#' `slope.weighting`, as in the table below. Furthermore, each subtype can
#' be either plain or smoothed (figure 4.44); smoothing is enabled setting the
#' logical argument `smoothed`.  (Not all combinations of arguments make
#' sense.)
#' 
#' ```
#'   Subtype | Rule       | Norm | Unbiased 
#'   --------|------------|------|---------
#'      a    | min step   |  --  |   NO 
#'      b    | max step   |  --  |   NO 
#'      c    | Di step    |   N  |  YES 
#'      d    | Di+Dj step | N+M  |  YES 
#' ```
#' 
#' 
#' **3. The Sakoe-Chiba set**
#' 
#' Sakoe-Chiba (Sakoe1978) discuss a family of slope-constrained patterns; they
#' are implemented as shown in page 47, table I.  Here, they are called
#' `symmetricP<x>` and `asymmetricP<x>`, where `<x>` corresponds
#' to Sakoe's integer slope parameter *P*.  Values available are
#' accordingly: `0` (no constraint), `1`, `05` (one half) and
#' `2`. See (Sakoe1978) for details.
#' 
#' **4. The Rabiner-Myers set**
#' 
#' The `type<XX><y>` step patterns follow the older Rabiner-Myers'
#' classification proposed in (Myers1980) and (MRR1980). Note that this is a
#' subset of the Rabiner-Juang set (Rabiner1993), and the latter should be
#' preferred in order to avoid confusion. `<XX>` is a Roman numeral
#' specifying the shape of the transitions; `<y>` is a letter in the range
#' `a-d` specifying the weighting used per step, as above; `typeIIx`
#' patterns also have a version ending in `s`, meaning the smoothing is
#' used (which does not permit skipping points). The `typeId, typeIId` and
#' `typeIIds` are unbiased and symmetric.
#' 
#' **5. Others**
#' 
#' The `rigid` pattern enforces a fixed unitary slope. It only makes sense
#' in combination with `open.begin=TRUE`, `open.end=TRUE` to find gapless
#' subsequences. It may be seen as the `P->inf` limiting case in Sakoe's classification.
#' 
#' `mori2006` is Mori's asymmetric step-constrained pattern (Mori2006). It
#' is normalized by the matched reference length.
#' 
#' [mvmStepPattern()] implements Latecki's Minimum Variance
#' Matching algorithm, and it is described in its own page.
#' 
#' 
#' **Methods**
#' 
#' `print.stepPattern` prints an user-readable description of the
#' recurrence equation defined by the given pattern.
#' 
#' `plot.stepPattern` graphically displays the step patterns productions
#' which can lead to element (0,0). Weights are shown along the step leading to
#' the corresponding element.
#' 
#' `t.stepPattern` transposes the productions and normalization hint so
#' that roles of query and reference become reversed.
#'
#'
#' @aliases stepPattern is.stepPattern print.stepPattern t.stepPattern
#' plot.stepPattern symmetric1 symmetric2 asymmetric rabinerJuangStepPattern
#' symmetricP0 asymmetricP0 symmetricP05 asymmetricP05 symmetricP1 asymmetricP1
#' symmetricP2 asymmetricP2 typeIa typeIas typeIb typeIbs typeIc typeIcs typeId
#' typeIds typeIIa typeIIb typeIIc typeIId typeIIIc typeIVc mori2006 rigid
#' @export symmetric1 symmetric2 asymmetric rabinerJuangStepPattern  symmetricP0 asymmetricP0 symmetricP05 asymmetricP05 symmetricP1 asymmetricP1  symmetricP2 asymmetricP2 typeIa typeIas typeIb typeIbs typeIc typeIcs typeId typeIds typeIIa typeIIb typeIIc typeIId typeIIIc typeIVc mori2006 rigid
#' @param x a step pattern object
#' @param type path specification, integer 1..7 (see (Rabiner1993), table 4.5)
#' @param slope.weighting slope weighting rule: character `"a"` to `"d"` (see (Rabiner1993), sec. 4.7.2.5)
#' @param smoothed logical, whether to use smoothing (see (Rabiner1993), fig. 4.44)
#' @param ... additional arguments to [print()].
#' @note Constructing `stepPattern` objects is tricky and thus undocumented. For a commented example please see source code for  `symmetricP1`.
#' @author Toni Giorgino
#' @seealso [mvmStepPattern()], implementing Latecki's Minimal
#' Variance Matching algorithm.
#' @references
#' * (GiorginoJSS) Toni Giorgino. *Computing and Visualizing
#' Dynamic Time Warping Alignments in R: The dtw Package.* Journal of
#' Statistical Software, 31(7), 1-24. \doi{10.18637/jss.v031.i07}
#' * (Itakura1975) Itakura, F., *Minimum prediction residual
#' principle applied to speech recognition,* Acoustics, Speech, and Signal
#' Processing, IEEE Transactions on , vol.23, no.1, pp.  67-72, Feb 1975. 
#' \doi{10.1109/TASSP.1975.1162641}
#' * (MRR1980) Myers,
#' C.; Rabiner, L. & Rosenberg, A. *Performance tradeoffs in dynamic time
#' warping algorithms for isolated word recognition*, IEEE Trans. Acoust.,
#' Speech, Signal Process., 1980, 28, 623-635. 
#' \doi{10.1109/TASSP.1980.1163491}
#' * (Mori2006) Mori,
#' A.; Uchida, S.; Kurazume, R.; Taniguchi, R.; Hasegawa, T. & Sakoe, H. Early
#' Recognition and Prediction of Gestures Proc. 18th International Conference
#' on Pattern Recognition ICPR 2006, 2006, 3, 560-563. 
#' \doi{10.1109/ICPR.2006.467}
#' * (Myers1980) Myers,
#' Cory S.  *A Comparative Study Of Several Dynamic Time Warping
#' Algorithms For Speech Recognition*, MS and BS thesis, Dept. of Electrical
#' Engineering and Computer Science, Massachusetts Institute of Technology,
#' archived Jun 20 1980, <https://hdl.handle.net/1721.1/27909> 
#' * (Rabiner1993) Rabiner, L. R., & Juang, B.-H. (1993). *Fundamentals of
#' speech recognition.* Englewood Cliffs, NJ: Prentice Hall. 
#' * (Sakoe1978) Sakoe, H.; Chiba, S., *Dynamic programming algorithm
#' optimization for spoken word recognition,* Acoustics, Speech, and Signal
#' Processing, IEEE Transactions on , vol.26, no.1, pp. 43-49, Feb 1978
#' \doi{10.1109/TASSP.1978.1163055}
#' @keywords ts
#' @examples
#' 
#' 
#' #########
#' ##
#' ## The usual (normalizable) symmetric step pattern
#' ## Step pattern recursion, defined as:
#' ## g[i,j] = min(
#' ##      g[i,j-1] + d[i,j] ,
#' ##      g[i-1,j-1] + 2 * d[i,j] ,
#' ##      g[i-1,j] + d[i,j] ,
#' ##   )
#' 
#' print(symmetric2)   # or just "symmetric2"
#' 
#' 
#' 
#' #########
#' ##
#' ## The well-known plotting style for step patterns
#' 
#' plot(symmetricP2,main="Sakoe's Symmetric P=2 recursion")
#' 
#' 
#' 
#' #########
#' ##
#' ## Same example seen in ?dtw , now with asymmetric step pattern
#' 
#' idx<-seq(0,6.28,len=100);
#' query<-sin(idx)+runif(100)/10;
#' reference<-cos(idx);
#' 
#' ## Do the computation 
#' asy<-dtw(query,reference,keep=TRUE,step=asymmetric);
#' 
#' dtwPlot(asy,type="density",main="Sine and cosine, asymmetric step")
#' 
#' 
#' #########
#' ##
#' ##  Hand-checkable example given in [Myers1980] p 61
#' ##
#' 
#' `tm` <-
#' structure(c(1, 3, 4, 4, 5, 2, 2, 3, 3, 4, 3, 1, 1, 1, 3, 4, 2,
#' 3, 3, 2, 5, 3, 4, 4, 1), .Dim = c(5L, 5L))
#' 
#' @name stepPattern
NULL

stepPattern <- function(v,norm=NA) {
  obj <- NULL;
  if(is.vector(v)) {
    obj <- matrix(v,ncol=4,byrow=TRUE);
  } else if(is.matrix(v)) {
    obj <- v;
  } else {
    stop("stepPattern constructor only supports vector or matrix");
  }
  class(obj)<-"stepPattern";
  attr(obj,"npat") <- max(obj[,1]);
  attr(obj,"norm") <- norm;
  return(obj);
}

#' @export
is.stepPattern <- function(x) {
  return(inherits(x,"stepPattern"));
}



## Transpose - exchange role of query and reference
#' @rdname stepPattern
#' @export
t.stepPattern <- function(x) {

  # exchange dx <-> dy
  tsp <- x[,c(1,3,2,4)];
  tsp <- stepPattern(tsp);

  # fix normalization, if available
  on <- attr(x,"norm");
  if(! is.na(on) ) {
    if(on == "N") {
      attr(tsp,"norm") <- "M";
    } else if(on == "M") {
      attr(tsp,"norm") <- "N";
    }
  }

  return(tsp);
}


## plot the step pattern
#' @rdname stepPattern
#' @export
plot.stepPattern <- function(x,...) {
  pats <- unique(x[,1]);                #list of patterns
  xr <- max(x[,2]);
  yr <- max(x[,3]);

  #for weight labels
  fudge <- c(-.5,1.2);                         
  alpha <- .5;                          # 1 start, 0 end

  ## dummy plot to fix the plot limits
  plot(-x[,2],-x[,3],type="n",
       xlab="Query index",ylab="Reference index",
       asp=1,lab=c(xr+1,yr+1,1),
       ax=FALSE,
       ...);

  for(i in pats) {
    ss <- x[,1]==i;
    lines(-x[ss,2],-x[ss,3],type="o", ...);

    if(sum(ss)==1) {
      next;                         
    }                               
    

    xh <- alpha*utils::head(x[ss,2],-1) + (1-alpha)*x[ss,2][-1];
    yh <- alpha*utils::head(x[ss,3],-1) + (1-alpha)*x[ss,3][-1];

    text(-xh,-yh,
         labels=round(x[ss,4][-1],2),
         adj=fudge,
         ...);
  }

  axis(1,at=c(-xr:0), ...)
  axis(2,at=c(-yr:0), ...)

  endpts <- x[,4]==-1;
  points(-x[endpts,2],-x[endpts,3],pch=16, ...);
}




## pretty-print the matrix meaning,
## so it will not be as write-only as now
#' @rdname stepPattern
#' @export
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
    steps<-matrix(steps[ns:1,],ncol=3); # enforce a matrix

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
        ccs<-ifelse(cc==1,"    ",sprintf("%2.2g *",cc));
        ds<-sprintf("+%s d[%s]",ccs,dijs);
        body<-paste(body,ds);
      }

    }

    body<-paste(body,",\n",s="");
  }

  tail<-")\n\n";

  norm <- attr(x,"norm");
  ntxt <- sprintf("Normalization hint: %s\n",norm);

  rv<-paste(head,body,tail,ntxt);

  cat("Step pattern recursion:\n");
  cat(rv);

}



## TODO: sanity check on the step pattern definition

.checkpattern <- function(sp) {
  ## must have 4 x n elements
  ## all integers
  ## first column in ascending order from 1, no missing steps
  ## 2nd, 3rd row non-negative
  ## 4th: first  for each step is -1
}


# Auxiliary function to easily map pattern -> delta

.mkDirDeltas <- function(dir) {
  m1 <- dir[ dir[,4]==-1, ,drop=FALSE ];
  m1 <- m1[,-4];
  m1 <- m1[,-1];
  return(m1);
}


## Extract rows belonging to pattern no. sn
## with first element stripped
## in reverse order

.extractpattern <- function(sp,sn) {
	sbs<-sp[,1]==sn;	# pick only rows beginning by sn
	spl<-sp[sbs,-1,drop=FALSE];
                                # of those: take only column Di, Dj, cost
                                # (drop first - pattern no. column)

	nr<-nrow(spl);	# how many are left
        spl<-spl[nr:1,,drop=FALSE];	# invert row order

	return(spl);
}





##################################################
##################################################


## Utility inner functions to manipulate
## step patterns. Could be implemented as
## a grammar, a'la ggplot2

.Pnew <- function(p,subt,smoo) {
  sp <- list();
  sp$i <- 0;
  sp$j <- 0;
  sp$p <- p;
  sp$subt <- subt;
  sp$smoo <- smoo;
  return(sp);
}

.Pstep <- function(sp,di,dj) {
  sp$i <- c(sp$i,di);
  sp$j <- c(sp$j,dj);
  return(sp);
}

.Pend <- function(sp,subt,smoo) {
  sp$si <- cumsum(sp$i);
  sp$sj <- cumsum(sp$j);
  sp$ni <- max(sp$si)-sp$si;
  sp$nj <- max(sp$sj)-sp$sj;

  w <- NULL;

  # smallest of i,j jumps
  if(sp$subt=="a") {
    w <- pmin(sp$i,sp$j);
  } else if(sp$subt=="b") {
    # largest of Di, Dj
    w <- pmax(sp$i,sp$j);
  } else if(sp$subt=="c") {
    # Di exactly
    w <- sp$i;
  } else if(sp$subt=="d") {
    # Di+Dj
    w <- sp$i+sp$j;
  } else {
    stop("Unsupported subtype");
  }

                                        # drop first element in w
  w <- w[-1];
  
  if(sp$smoo)
    w <- rep(mean(w),length(w));

                                        # prepend -1
  w <- c(-1,w);
  sp$w <- w;

  return(sp);
}

.PtoMx <- function(sp) {
  nr <- length(sp$i);
  mx <- matrix(nrow=nr,ncol=4)
  mx[,1] <- sp$p;
  mx[,2] <- sp$ni;
  mx[,3] <- sp$nj;
  mx[,4] <- sp$w;
  return(mx);
}


#' @export
#' @rdname stepPattern
rabinerJuangStepPattern <- function(type,slope.weighting="d",smoothed=FALSE) {

  sw <- slope.weighting;
  sm <- smoothed;
  
  ## Actually build the step
  r <- switch(type,
              .RJtypeI(sw,sm),
              .RJtypeII(sw,sm),
              .RJtypeIII(sw,sm),
              .RJtypeIV(sw,sm),
              .RJtypeV(sw,sm),
              .RJtypeVI(sw,sm),
              .RJtypeVII(sw,sm)
              );

  norm <- NA;
  if(sw=="c") {
    norm <- "N";
  } else if(sw=="d") {
    norm <- "N+M";
  }

  # brain-damaged legacy
  rv <- as.vector(t(r));                
  rs <- stepPattern(rv);
  attr(rs,"norm") <- norm;
  attr(rs,"call") <- match.call();
  return(rs);
}



.RJtypeI <- function(s,m) {
  t <- .Pnew(1,s,m)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m1 <- .PtoMx(t);

  t <- .Pnew(2,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pend(t);
  m2 <- .PtoMx(t);

  t <- .Pnew(3,s,m)
  t <- .Pstep(t,0,1)
  t <- .Pend(t);
  m3 <- .PtoMx(t)

  return(rbind(m1,m2,m3));
}


.RJtypeII <- function(s,m) {
  t <- .Pnew(1,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m1 <- .PtoMx(t);

  t <- .Pnew(2,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pend(t);
  m2 <- .PtoMx(t);

  t <- .Pnew(3,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,0,1)
  t <- .Pend(t);
  m3 <- .PtoMx(t)

  return(rbind(m1,m2,m3));
}



.RJtypeIII <- function(s,m) {
  t <- .Pnew(1,s,m)
  t <- .Pstep(t,2,1)
  t <- .Pend(t);
  m1 <- .PtoMx(t);

  t <- .Pnew(2,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pend(t);
  m2 <- .PtoMx(t);

  t <- .Pnew(3,s,m)
  t <- .Pstep(t,1,2)
  t <- .Pend(t);
  m3 <- .PtoMx(t)

  return(rbind(m1,m2,m3));
}



.RJtypeIV <- function(s,m) {
  t <- .Pnew(1,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m1 <- .PtoMx(t);

  t <- .Pnew(2,s,m)
  t <- .Pstep(t,1,2)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m2 <- .PtoMx(t);

  t <- .Pnew(3,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pend(t);
  m3 <- .PtoMx(t)

  t <- .Pnew(4,s,m)
  t <- .Pstep(t,1,2)
  t <- .Pend(t);
  m4 <- .PtoMx(t)

  return(rbind(m1,m2,m3,m4));
}



.RJtypeV <- function(s,m) {
  t <- .Pnew(1,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,0)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m1 <- .PtoMx(t);

  t <- .Pnew(2,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m2 <- .PtoMx(t);

  t <- .Pnew(3,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pend(t);
  m3 <- .PtoMx(t)

  t <- .Pnew(4,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,0,1)
  t <- .Pend(t);
  m4 <- .PtoMx(t)

  t <- .Pnew(5,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,0,1)
  t <- .Pstep(t,0,1)
  t <- .Pend(t);
  m5 <- .PtoMx(t)

  return(rbind(m1,m2,m3,m4,m5));
}



.RJtypeVI <- function(s,m) {
  t <- .Pnew(1,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m1 <- .PtoMx(t);

  t <- .Pnew(2,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pend(t);
  m2 <- .PtoMx(t);

  t <- .Pnew(3,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,0,1)
  t <- .Pend(t);
  m3 <- .PtoMx(t)

  return(rbind(m1,m2,m3));
}




.RJtypeVII <- function(s,m) {
  t <- .Pnew(1,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,0)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m1 <- .PtoMx(t);

  t <- .Pnew(2,s,m)
  t <- .Pstep(t,1,2)
  t <- .Pstep(t,1,0)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m2 <- .PtoMx(t);

  t <- .Pnew(3,s,m)
  t <- .Pstep(t,1,3)
  t <- .Pstep(t,1,0)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m3 <- .PtoMx(t)

  t <- .Pnew(4,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m4 <- .PtoMx(t)

  t <- .Pnew(5,s,m)
  t <- .Pstep(t,1,2)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m5 <- .PtoMx(t)

  t <- .Pnew(6,s,m)
  t <- .Pstep(t,1,3)
  t <- .Pstep(t,1,0)
  t <- .Pend(t);
  m6 <- .PtoMx(t);

  t <- .Pnew(7,s,m)
  t <- .Pstep(t,1,1)
  t <- .Pend(t);
  m7 <- .PtoMx(t)

  t <- .Pnew(8,s,m)
  t <- .Pstep(t,1,2)
  t <- .Pend(t);
  m8 <- .PtoMx(t)

  t <- .Pnew(9,s,m)
  t <- .Pstep(t,1,3)
  t <- .Pend(t);
  m9 <- .PtoMx(t)

  return(rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9));
}














##################################################
##################################################


##
## Various step patterns, defined as internal variables
##
## First column: enumerates step patterns.
## Second   	 step in query index
## Third	 step in reference index
## Fourth	 weight if positive, or -1 if starting point
##
## For \cite{} see dtw.bib in the package
##



## Widely-known variants

## White-Neely symmetric (default)
## aka Quasi-symmetric \cite{White1976}
## normalization: no (N+M?)
symmetric1 <- stepPattern(c(
                            1,1,1,-1,
                            1,0,0,1,
                            2,0,1,-1,
                            2,0,0,1,
                            3,1,0,-1,
                            3,0,0,1
                            ));


## Normal symmetric
## normalization: N+M
symmetric2 <- stepPattern(c(
                            1,1,1,-1,
                            1,0,0,2,
                            2,0,1,-1,
                            2,0,0,1,
                            3,1,0,-1,
                            3,0,0,1
                            ),"N+M");


## classic asymmetric pattern: max slope 2, min slope 0
## normalization: N
asymmetric <-  stepPattern(c(
                             1,1,0,-1,
                             1,0,0,1,
                             2,1,1,-1,
                             2,0,0,1,
                             3,1,2,-1,
                             3,0,0,1
                           ),"N");


# % \item{\code{symmetricVelichkoZagoruyko}}{symmetric, reproduced from %
# [Sakoe1978]. Use distance matrix \code{1-d}}
# 

## normalization: max[N,M]
## note: local distance matrix is 1-d
## \cite{Velichko}
.symmetricVelichkoZagoruyko <- stepPattern(c(
		1, 0, 1, -1,
		2, 1, 1, -1,
		2, 0, 0, -1.001,
		3, 1, 0, -1 ));


# % \item{\code{asymmetricItakura}}{asymmetric, slope contrained 0.5 -- 2
# from reference [Itakura1975]. This is the recursive definition % that
# generates the Itakura parallelogram; }
# 

## Itakura slope-limited asymmetric \cite{Itakura1975}
## Max slope: 2; min slope: 1/2
## normalization: N
.asymmetricItakura <-  stepPattern(c(
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







#############################
## Slope-limited versions
##
## Taken from Table I, page 47 of "Dynamic programming algorithm
## optimization for spoken word recognition," Acoustics, Speech, and
## Signal Processing, vol.26, no.1, pp. 43-49, Feb 1978 URL:
## http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055
##
## Mostly unchecked



## Row P=0
symmetricP0 <- symmetric2;

## normalization: N ?
asymmetricP0 <- stepPattern(c(
                                  1,0,1,-1,
                                  1,0,0,0,
                                  2,1,1,-1,
                                  2,0,0,1,
                                  3,1,0,-1,
                                  3,0,0,1
                                ),"N");


## alternative implementation
.asymmetricP0b <- stepPattern(c(
                                  1,0,1,-1,
                                  2,1,1,-1,
                                  2,0,0,1,
                                  3,1,0,-1,
                                  3,0,0,1
                                ),"N");



## Row P=1/2
symmetricP05 <-  stepPattern(c(
                        1  ,  1, 3 , -1,
                        1  ,  0, 2 ,  2,
                        1  ,  0, 1 ,  1,
                        1  ,  0, 0 ,  1,
                        2  ,  1, 2 , -1,
                        2  ,  0, 1 ,  2,
                        2  ,  0, 0 ,  1,
                        3  ,  1, 1 , -1,
                        3  ,  0, 0 ,  2,
                        4  ,  2, 1 , -1,
                        4  ,  1, 0 ,  2,
                        4  ,  0, 0 ,  1,
                        5  ,  3, 1 , -1,
                        5  ,  2, 0 ,  2,
                        5  ,  1, 0 ,  1,
                        5  ,  0, 0 ,  1
                               ),"N+M");

asymmetricP05 <-  stepPattern(c(
                        1  , 1 , 3 , -1,
                        1  , 0 , 2 ,1/3,
                        1  , 0 , 1 ,1/3,
                        1  , 0 , 0 ,1/3,
                        2  , 1 , 2 , -1,
                        2  , 0 , 1 , .5,
                        2  , 0 , 0 , .5,
                        3  , 1 , 1 , -1,
                        3  , 0 , 0 , 1 ,
                        4  , 2 , 1 , -1,
                        4  , 1 , 0 , 1 ,
                        4  , 0 , 0 , 1 ,
                        5  , 3 , 1 , -1,
                        5  , 2 , 0 , 1 ,
                        5  , 1 , 0 , 1 ,
                        5  , 0 , 0 , 1
                               ),"N");



## Row P=1
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
                        ),"N+M");

asymmetricP1 <- stepPattern(c(
                              1, 1 , 2 , -1 ,
                              1, 0 , 1 , .5 ,
                              1, 0 , 0 , .5 ,
                              2, 1 , 1 , -1 ,
                              2, 0 , 0 ,  1 ,
                              3, 2 , 1 , -1 ,
                              3, 1 , 0 ,  1 ,
                              3, 0 , 0 ,  1
                              ),"N");


## Row P=2
symmetricP2 <- stepPattern(c(
	1, 2, 3, -1,
	1, 1, 2, 2,
	1, 0, 1, 2,
	1, 0, 0, 1,
	2, 1, 1, -1,
	2, 0, 0, 2,
	3, 3, 2, -1,
	3, 2, 1, 2,
	3, 1, 0, 2,
	3, 0, 0, 1
),"N+M");

asymmetricP2 <- stepPattern(c(
	1, 2 , 3  , -1,
	1, 1 , 2  ,2/3,
	1, 0 , 1  ,2/3,
	1, 0 , 0  ,2/3,
	2, 1 , 1  ,-1 ,
	2, 0 , 0  ,1  ,
	3, 3 , 2  ,-1 ,
	3, 2 , 1  ,1  ,
	3, 1 , 0  ,1  ,
	3, 0 , 0  ,1
),"N");






################################
## Taken from Table III, page 49.
## Four varieties of DP-algorithm compared

## 1st row:  asymmetric

## 2nd row:  symmetricVelichkoZagoruyko

## 3rd row:  symmetric1

## 4th row:  asymmetricItakura




#############################
## Classified according to Rabiner
##
## Taken from chapter 2, Myers' thesis [4]. Letter is
## the weighting function:
##
##      rule       norm   unbiased
##   a  min step   ~N     NO
##   b  max step   ~N     NO
##   c  x step     N      YES
##   d  x+y step   N+M    YES
##
## Mostly unchecked

# R-Myers     R-Juang
# type I      type II   
# type II     type III
# type III    type IV
# type IV     type VII


typeIa <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0,  1,
                         1, 0, 0,  0,
                         2, 1, 1, -1,
                         2, 0, 0,  1,
                         3, 1, 2, -1,
                         3, 0, 1,  1,
                         3, 0, 0,  0
 ));

typeIb <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0,  1,
                         1, 0, 0,  1,
                         2, 1, 1, -1,
                         2, 0, 0,  1,
                         3, 1, 2, -1,
                         3, 0, 1,  1,
                         3, 0, 0,  1
 ));

typeIc <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0,  1,
                         1, 0, 0,  1,
                         2, 1, 1, -1,
                         2, 0, 0,  1,
                         3, 1, 2, -1,
                         3, 0, 1,  1,
                         3, 0, 0,  0
 ),"N");

typeId <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0,  2,
                         1, 0, 0,  1,
                         2, 1, 1, -1,
                         2, 0, 0,  2,
                         3, 1, 2, -1,
                         3, 0, 1,  2,
                         3, 0, 0,  1
 ),"N+M");

## ----------
## smoothed variants of above

typeIas <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0, .5,
                         1, 0, 0, .5,
                         2, 1, 1, -1,
                         2, 0, 0,  1,
                         3, 1, 2, -1,
                         3, 0, 1, .5,
                         3, 0, 0, .5
 ));


typeIbs <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0,  1,
                         1, 0, 0,  1,
                         2, 1, 1, -1,
                         2, 0, 0,  1,
                         3, 1, 2, -1,
                         3, 0, 1,  1,
                         3, 0, 0,  1
 ));


typeIcs <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0,  1,
                         1, 0, 0,  1,
                         2, 1, 1, -1,
                         2, 0, 0,  1,
                         3, 1, 2, -1,
                         3, 0, 1, .5,
                         3, 0, 0, .5
 ),"N");


typeIds <-  stepPattern(c(
                         1, 2, 1, -1,
                         1, 1, 0, 1.5,
                         1, 0, 0, 1.5,
                         2, 1, 1, -1,
                         2, 0, 0,  2,
                         3, 1, 2, -1,
                         3, 0, 1, 1.5,
                         3, 0, 0, 1.5
 ),"N+M");






## ----------

typeIIa <- stepPattern(c(
                        1,  1,  1, -1,
                        1,  0,  0, 1,
                        2,  1,  2, -1,
                        2,  0,  0, 1,
                        3,  2,  1, -1,
                        3,  0,  0, 1
                        ));

typeIIb <- stepPattern(c(
                        1,  1,  1, -1,
                        1,  0,  0, 1,
                        2,  1,  2, -1,
                        2,  0,  0, 2,
                        3,  2,  1, -1,
                        3,  0,  0, 2
                        ));

typeIIc <- stepPattern(c(
                        1,  1,  1, -1,
                        1,  0,  0, 1,
                        2,  1,  2, -1,
                        2,  0,  0, 1,
                        3,  2,  1, -1,
                        3,  0,  0, 2
                        ),"N");

typeIId <- stepPattern(c(
                        1,  1,  1, -1,
                        1,  0,  0, 2,
                        2,  1,  2, -1,
                        2,  0,  0, 3,
                        3,  2,  1, -1,
                        3,  0,  0, 3
                        ),"N+M");

## ----------

## Rabiner [3] discusses why this is not equivalent to Itakura's

typeIIIc <-  stepPattern(c(
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
                       ),"N");



## ----------

## numbers follow as production rules in fig 2.16

typeIVc <-  stepPattern(c(
                          1,  1,  1,  -1,
                          1,  0,  0,   1,
                          2,  1,  2,  -1,
                          2,  0,  0,   1,
                          3,  1,  3,  -1,
                          3,  0,  0,   1,
                          4,  2,  1,  -1,
                          4,  1,  0,   1,
                          4,  0,  0,   1,
                          5,  2,  2,  -1,
                          5,  1,  0,   1,
                          5,  0,  0,   1,
                          6,  2,  3,  -1,
                          6,  1,  0,   1,
                          6,  0,  0,   1,
                          7,  3,  1,  -1,
                          7,  2,  0,   1,
                          7,  1,  0,   1,
                          7,  0,  0,   1,
                          8,  3,  2,  -1,
                          8,  2,  0,   1,
                          8,  1,  0,   1,
                          8,  0,  0,   1,
                          9,  3,  3,  -1,
                          9,  2,  0,   1,
                          9,  1,  0,   1,
                          9,  0,  0,   1
 ),"N");






#############################
## 
## Mori's asymmetric step-constrained pattern. Normalized in the
## reference length.
##
## Mori, A.; Uchida, S.; Kurazume, R.; Taniguchi, R.; Hasegawa, T. &
## Sakoe, H. Early Recognition and Prediction of Gestures Proc. 18th
## International Conference on Pattern Recognition ICPR 2006, 2006, 3,
## 560-563
##

mori2006 <-  stepPattern(c(
                           1, 2, 1, -1,
                           1, 1, 0,  2,
                           1, 0, 0,  1,
                           2, 1, 1, -1,
                           2, 0, 0,  3,
                           3, 1, 2, -1,
                           3, 0, 1,  3,
                           3, 0, 0,  3
 ),"M");


## Completely unflexible: fixed slope 1. Only makes sense with
## open.begin and open.end
rigid <- stepPattern(c(1,1,1,-1,
                       1,0,0,1  ),"N")


\name{stepPattern}

\alias{stepPattern}
\alias{is.stepPattern}
\alias{print.stepPattern}
\alias{t.stepPattern}
\alias{plot.stepPattern}

\alias{symmetric1}
\alias{symmetric2}
\alias{asymmetric}

%\alias{asymmetricItakura}
%\alias{symmetricVelichkoZagoruyko}

\alias{rabinerJuangStepPattern}

\alias{symmetricP0}
\alias{asymmetricP0}
\alias{symmetricP05}
\alias{asymmetricP05}
\alias{symmetricP1}
\alias{asymmetricP1}
\alias{symmetricP2}
\alias{asymmetricP2}


\alias{typeIa}
\alias{typeIas}
\alias{typeIb}
\alias{typeIbs} 
\alias{typeIc}
\alias{typeIcs}
\alias{typeId}
\alias{typeIds}
\alias{typeIIa}
\alias{typeIIb}
\alias{typeIIc}
\alias{typeIId}
\alias{typeIIIc} 
\alias{typeIVc}

\alias{mori2006}



\title{Step patterns for DTW}

\description{ A \code{stepPattern} object lists the transitions
 allowed while searching for the minimum-distance path.  DTW variants
 are implemented by passing one of the objects described in this page
 to the \code{stepPattern} argument of the \code{\link{dtw}}  call.  }

\usage{
## Well-known step patterns
  symmetric1
  symmetric2
  asymmetric
%  asymmetricItakura
%  symmetricVelichkoZagoruyko


## Step patterns classified according to Rabiner-Juang [3]
  rabinerJuangStepPattern(type,slope.weighting="d",smoothed=FALSE)

## Slope-constrained step patterns from Sakoe-Chiba [1]
  symmetricP0;  asymmetricP0
  symmetricP05; asymmetricP05
  symmetricP1;  asymmetricP1
  symmetricP2;  asymmetricP2

## Step patterns classified according to Rabiner-Myers [4]
  typeIa;   typeIb;   typeIc;   typeId;
  typeIas;  typeIbs;  typeIcs;  typeIds;  # smoothed
  typeIIa;  typeIIb;  typeIIc;  typeIId;
  typeIIIc; typeIVc;

## Miscellaneous
  mori2006;



\method{print}{stepPattern}(x,...)
\method{plot}{stepPattern}(x,...)
\method{t}{stepPattern}(x)

stepPattern(v,norm=NA)
is.stepPattern(x)

}

\arguments{
  \item{x}{a step pattern object}
  \item{type}{path specification, integer 1..7 (see [3], table 4.5)}
  \item{slope.weighting}{slope weighting rule: character \code{"a"}
    to \code{"d"} (see [3],    sec. 4.7.2.5)}
  \item{smoothed}{logical, whether to use smoothing (see [3],
  fig. 4.44) }
  \item{v}{a vector defining the stepPattern structure}
  \item{norm}{normalization hint (character)}
  \item{...}{additional arguments to \code{\link{print}}.}
}


\details{

  A step pattern characterizes the matching model and slope constraint
  specific of a DTW variant. They also known as local- or
  slope-constraints, transition types, production or recursion rules
  [7].

  \code{print.stepPattern} prints an user-readable
  description of the recurrence equation defined by the given pattern.

  \code{plot.stepPattern} graphically displays the step patterns
  productions which can lead to element (0,0). Weights are 
  shown along the step leading to the corresponding element.

  \code{t.stepPattern} transposes the productions and normalization hint
  so that roles of query and reference become reversed.
  

  A variety of classifications have been proposed for step patterns,
  including Sakoe-Chiba [1]; Rabiner-Juang [3]; and Rabiner-Myers [4].
  The \code{dtw} package implements all of the transition types found in
  those papers, with the exception of Itakura's and Velichko-Zagoruyko's
  steps which require subtly different algorithms (this may be rectified
  in the future). Itakura recursion is almost, but not quite, equivalent
  to \code{typeIIIc}.
  
  For convenience, we shall review pre-defined step patterns grouped by
  classification. Note that the same pattern may be listed under
  different names. Refer to paper [7] for full details.


  
\strong{1. Well-known step patterns}

  These common transition types are used in quite a lot of implementations.

  \code{symmetric1} (or White-Neely) is the commonly used
  quasi-symmetric, no local constraint, non-normalizable. It is biased
  in favor of oblique steps.
  
  \code{symmetric2} is normalizable, symmetric, with no local slope
  constraints.  Since one diagonal step costs as much as the two
  equivalent steps along the sides, it can be normalized dividing by
  \code{N+M} (query+reference lengths).

  \code{asymmetric} is asymmetric, slope constrained between 0 and
  2. Matches each element of the query time series exactly once, so
  the warping path \code{index2~index1} is guaranteed to
  be single-valued.  Normalized by \code{N} (length of query).

%    \item{\code{asymmetricItakura}}{asymmetric, slope contrained 0.5
%	-- 2 from reference [2]. This is the recursive definition
%	that generates the Itakura parallelogram; }

%    \item{\code{symmetricVelichkoZagoruyko}}{symmetric, reproduced from
%    [1]. Use distance matrix \code{1-d}}



\strong{2. The Rabiner-Juang set}

  A comprehensive table of step patterns is proposed by Rabiner-Juang
  [3], tab. 4.5.  All of them can be recovered by the
  \code{rabinerJuangStepPattern(type,slope.weighting,smoothed)}
  function.

  Seven families, labelled with Roman numerals I-VII, are
  selected through the integer argument \code{type}. Each family has
  four slope weighting sub-types, named in sec. 4.7.2.5 as "Type (a)" to
  "Type (d)"; they are selected passing a character argument
  \code{slope.weighting}, as in the table below. Furthermore, each
  subtype can be plain or smoothed (figure 4.44); smoothing is enabled
  setting the logical argument \code{smoothed}.  (Not all combinations
  of arguments make sense.)

  \tabular{cccc}{
     Subtype \tab    Rule  \tab     Norm  \tab   Unbiased    \cr
%    --------------------------------
        a    \tab  min step \tab    --    \tab   NO          \cr
        b    \tab  max step \tab    --    \tab   NO          \cr
        c    \tab  Di step   \tab    N     \tab   YES         \cr
        d    \tab  Di+Dj step \tab    N+M   \tab   YES         \cr
    }
  



\strong{3. The Sakoe-Chiba set}

  \code{symmetricPx} is the family of Sakoe's symmetric steps, slope
  contraint \code{x}; \code{asymmetricPx} are Sakoe's asymmetric, slope
  contraint \code{x}. These slope-constrained patterns are discussed in
  Sakoe-Chiba [1], and implemented as shown in page 47, table I. Values
  available for \emph{P} (\code{x}) are accordingly: \code{0} (no
  constraint), \code{1}, \code{05} (one half) and \code{2}. See
  reference for details.


  

\strong{4. The Rabiner-Myers set}
  
  The \code{typeXXx} step patterns follow the older Rabiner-Myers'
  classification given in [4-5]. Note that they are a subset of the
  Rabiner-Juang set [3], which should be preferred to avoid
  confusion. \code{XX} is a roman numeral specifying the shape of the
  transitions; \code{x} is a letter in the range \code{a-d} according
  the type of weighting used per step, as above; \code{typeIIx} patterns
  also have a version ending in \code{s} meaning the path smoothing is
  used (which does not permit skipping points). The \code{typeId,
  typeIId} and \code{typeIIds} are unbiased and symmetric.



\strong{5. Other}

  Mori's [6] asymmetric step-constrained pattern is called
  \code{mori2006}. It is normalized in the reference length.


}


\note{    
  The \code{stepPattern} constructor is currently not well
  documented. For a commented example please see source code for
  \code{symmetricP1}.
}



\references{ 
  [1] Sakoe, H.; Chiba, S., \emph{Dynamic programming algorithm
  optimization for spoken word recognition,} Acoustics, Speech, and
  Signal Processing [see also IEEE Transactions on Signal Processing],
  IEEE Transactions on , vol.26, no.1, pp. 43-49, Feb 1978 URL:
  \url{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055}

  [2] Itakura, F., \emph{Minimum prediction residual principle applied
  to speech recognition,} Acoustics, Speech, and Signal Processing [see
  also IEEE Transactions on Signal Processing], IEEE Transactions on ,
  vol.23, no.1, pp.  67-72, Feb 1975. URL:
  \url{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1162641}

  [3] Rabiner, L. R., & Juang, B.-H. (1993). \emph{Fundamentals of speech
  recognition.} Englewood Cliffs, NJ: Prentice Hall.

  [4] Myers, C. S.  \emph{A Comparative Study Of Several Dynamic Time
  Warping Algorithms For Speech Recognition}, MS and BS thesis, MIT Jun
  20 1980, \url{dspace.mit.edu/bitstream/1721.1/27909/1/07888629.pdf}

  [5] Myers, C.; Rabiner, L. & Rosenberg, A. \emph{Performance tradeoffs in
  dynamic time warping algorithms for isolated word recognition},
  IEEE Trans. Acoust., Speech, Signal Process., 1980, 28, 623-635

  [6] Mori, A.; Uchida, S.; Kurazume, R.; Taniguchi, R.; Hasegawa, T. &
  Sakoe, H. Early Recognition and Prediction of Gestures Proc. 18th
  International Conference on Pattern Recognition ICPR 2006, 2006, 3,
  560-563
  
  [7] Toni Giorgino. \emph{Computing and Visualizing Dynamic Time Warping
  Alignments in R: The dtw Package.}  Journal of Statistical
  Software, 31(7), 1-24. \url{http://www.jstatsoft.org/v31/i07/}
}




\author{Toni Giorgino}

\examples{

#########
##
## The usual (normalizable) symmetric step pattern
## Step pattern recursion, defined as:
## g[i,j] = min(
##      g[i,j-1] + d[i,j] ,
##      g[i-1,j-1] + 2 * d[i,j] ,
##      g[i-1,j] + d[i,j] ,
##   )

print(symmetric2)   # or just "symmetric2"



#########
##
## The well-known plotting style for step patterns

plot(symmetricP2,main="Sakoe's Symmetric P=2 recursion")



#########
##
## Same example seen in ?dtw , now with asymmetric step pattern

idx<-seq(0,6.28,len=100);
query<-sin(idx)+runif(100)/10;
reference<-cos(idx);

## Do the computation 
asy<-dtw(query,reference,keep=TRUE,step=asymmetric);

dtwPlot(asy,type="density",main="Sine and cosine, asymmetric step")


#########
##
##  Hand-checkable example given in [4] p 61
##

`tm` <-
structure(c(1, 3, 4, 4, 5, 2, 2, 3, 3, 4, 3, 1, 1, 1, 3, 4, 2,
3, 3, 2, 5, 3, 4, 4, 1), .Dim = c(5L, 5L))



}



\concept{Dynamic Time Warp}
\concept{Dynamic Programming}
\concept{Step pattern}
\concept{Transition}
\concept{Local constraint}
\concept{Asymmetric DTW}
\concept{Symmetric DTW}



\keyword{ ts }

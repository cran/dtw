\name{stepPattern}

\alias{stepPattern}
\alias{is.stepPattern}
\alias{print.stepPattern}

\alias{symmetric1}
\alias{symmetric2}
\alias{asymmetric}

\alias{asymmetricItakura}
%\alias{symmetricVelichkoZagoruyko}

\alias{symmetricP0}
\alias{asymmetricP0}
\alias{symmetricP05}
\alias{asymmetricP05}
\alias{symmetricP1}
\alias{asymmetricP1}
\alias{symmetricP2}
\alias{asymmetricP2}




\title{Local constraints and step patterns for DTW}

\description{ DTW variants are implemented through step pattern objects.
  A \code{stepPattern} instance lists the transitions allowed by the
  \code{\link{dtw}} function in the search for the minimum-distance
  path.  }

\usage{
## Well-known step patterns
  symmetric1
  symmetric2
  asymmetric
  asymmetricItakura
%  symmetricVelichkoZagoruyko

## Slope-constrained step patterns from Sakoe-Chiba
  symmetricP0;  asymmetricP0
  symmetricP05; asymmetricP05
  symmetricP1;  asymmetricP1
  symmetricP2;  asymmetricP2 


\method{print}{stepPattern}(x,...)

stepPattern(v)
is.stepPattern(x)

}

\arguments{
  \item{x}{a step pattern object}
  \item{v}{a vector defining the stepPattern structure (see below)}
  \item{...}{additional arguments to \code{\link{print}}.}
}


\details{

  A step pattern characterizes the matching model and/or slope constraint
  specific of a DTW variant.

  \code{print.stepPattern} prints an user-readable
  description of the recurrence equation defined by the given pattern.


  % TODO REPLACE as TABLE
  Several step patterns are pre-defined with the package:

  \itemize{

    \item{\code{symmetric1}}{quasi-symmetric, no slope constraint
    (favours oblique steps);}

    \item{\code{symmetric2}}{properly symmetric, no slope constraint:
    one diagonal step costs as much as the sum of the equivalent steps
    along the sides; }

    \item{\code{asymmetric}}{asymmetric, slope constrained between 0 and
      2. Matches each element of the query time series exactly once, so
      the warping path \code{index2~index1} is guaranteed to
      be single-valued; }

    \item{\code{asymmetricItakura}}{asymmetric, slope contrained 0.5
	-- 2 from reference [2]. This is the recursive definition
	that generates the Itakura parallelogram; }

%    \item{\code{symmetricVelichkoZagoruyko}}{symmetric, reproduced from [1]. Use distance matrix \code{1-d}}

    \item{\code{symmetricPx}}{symmetric, slope contraint $P=x$;}

    \item{\code{asymmetricPx}}{asymmetric, slope contraint $P=x$.}


  }

  The \code{symmetricPx} and \code{asymmetricPx} slope-constrained
  patterns are discussed in Sakoe-Chiba [1], and reproduced as shown in
  page 47, table I. Values available for \emph{P} are accordingly: \code{0}
  (no constraint), \code{1}, \code{05} (half) and \code{2}. See [1] for
  details. They are also known as type V (Rabiner and Huang) [3].
  
  The \code{stepPattern} constructor is currently not well
  documented. Please see the example below, implementing Sakoe's
  \emph{P=1, Symmetric} algorithm.

    \preformatted{
     symmetricP1 <- stepPattern(c(
        1,1,2,-1,       # First branch: g(i-1,j-2) +
        1,0,1,2,        #            + 2d( i ,j-1) + 
        1,0,0,1,        #            +  d( i , j )
        2,1,1,-1,       # Second br.:   g(i-1,j-1) +
        2,0,0,2,        #            + 2d( i , j )
        3,2,1,-1,       # Third branch: g(i-2,j-1) +
        3,1,0,2,        #            + 2d(i-1, j ) +
        3,0,0,1         #            +  d( i , j )              
      ));
    }

    Decoding is left to the reader as an exercise, and
    \code{print.stepPattern} may come handy.
		

}



\references{ 
[1] Sakoe, H.; Chiba, S., \emph{Dynamic programming algorithm optimization for spoken word recognition,}
 Acoustics, Speech, and Signal Processing [see also IEEE Transactions on Signal Processing], IEEE Transactions on , 
vol.26, no.1, pp. 43-49, Feb 1978 URL: \url{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055} \cr
\cr
[2] Itakura, F., \emph{Minimum prediction residual principle applied to speech
recognition,} Acoustics, Speech, and Signal Processing [see also IEEE
Transactions on Signal Processing], IEEE Transactions on , vol.23, no.1, pp.
67-72, Feb 1975. URL:
\url{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1162641} \cr
\cr
[3] Rabiner, L. R., & Juang, B.-H. (1993). \emph{Fundamentals of speech
recognition.} Englewood Cliffs, NJ: Prentice Hall.
 
}




\author{Toni Giorgino}

\examples{

## The usual (normalizable) symmetric step pattern
## Step pattern recursion, defined as:
## g[i,j] = min(
##      g[i,j-1] + d[i,j] ,
##      g[i-1,j-1] + 2 * d[i,j] ,
##      g[i-1,j] + d[i,j] ,
##   )

print.stepPattern(symmetric2)   # or just "symmetric2"



## Same example seen in ?dtw , now with asymmetric step pattern

idx<-seq(0,6.28,len=100);
query<-sin(idx)+runif(100)/10;
template<-cos(idx);

## Do the computation 
asy<-dtw(query,template,keep=TRUE,step=asymmetric);

dtwPlot(asy,type="density",main="Sine and cosine, asymmetric step")

}



\concept{Dynamic Time Warp}
\concept{Dynamic Programming}
\concept{Step pattern}
\concept{Transition}
\concept{Local constraint}
\concept{Asymmetric DTW}
\concept{Symmetric DTW}



\keyword{ ts }
    
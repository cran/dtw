\name{stepPattern}

\alias{stepPattern}
\alias{is.stepPattern}
\alias{print.stepPattern}

\alias{symmetric1}
\alias{symmetric2}
\alias{asymmetric}
\alias{symmetricP1}
\alias{asymmetricItakura}


\title{Local constraints and step patterns for DTW}

\description{ A step pattern object lists the allowed step patterns for
a given DTW variant.  It can be used to constrain the warping
  paths considered by the   \code{\link{dtw}} function.
}

\usage{

# predefined step patterns:
symmetric1
symmetric2
asymmetric
symmetricP1
asymmetricItakura

\method{print}{stepPattern}(x,...)
is.stepPattern(x)
stepPattern(v)

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
    description of the recurrence equation defined by the pattern.

    Several step patterns are pre-defined:

    \itemize {
    \item{\code{symmetric1}}{quasi-symmetric, no slope constraint. Favours oblique steps.}

    \item{\code{symmetric2}}{properly symmetric, no slope constraint: oblique
    step costs as much as the sum of the equivalent steps along the
    sides.}

    \item{\code{asymmetric}}{asymmetric, no slope constraint. Matches each
      element of the query time series exactly once.}

    \item{\code{symmetricP1}}{symmetric, slope contraint $P=1$ (see
      reference [1], page 47, table I)}

    \item{\code{asymmetricItakura}}{asymmetric, slope contrained. This is the 
	recursive function described  in reference [2] that generates 
  	the Itakura parallelogram. 
    }

    \item{}{...more to come}

    }
    

    The vector description for \code{stepPattern} constructor is
    currently not well documented. Please see the example below,
    implementing Sakoe's P=1, Symmetric algorithm in Table I, p47 [1].

    \preformatted{
      symmetricP1 <- stepPattern(c(
                         1,1,2,-1,	# First branch: g(i-1,j-2) +
                         1,0,1,2,	#            + 2d( i ,j-1) + 
                         1,0,0,1,	#            +  d( i , j )
                         2,1,1,-1,	# Second branch: g(i-1,j-1) +
                         2,0,0,2,	#             + 2d( i , j )
                         3,2,1,-1,	# Third branch: g(i-2,j-1) +
                         3,1,0,2,	#            + 2d(i-1, j ) +
                         3,0,0,1	#            +  d( i , j )		
                         ));
			 }

			 Decoding is left
    to the reader as an exercise.  \code{print.stepPattern} may come handy.
		

}



\references{ 
[1] Sakoe, H.; Chiba, S., \emph{Dynamic programming algorithm optimization for spoken word recognition,}
 Acoustics, Speech, and Signal Processing [see also IEEE Transactions on Signal Processing], IEEE Transactions on , 
vol.26, no.1, pp. 43-49, Feb 1978 URL: \url{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055} \cr
\cr
[2] Itakura, F., \emph{Minimum prediction residual principle applied to speech
recognition,} Acoustics, Speech, and Signal Processing [see also IEEE
Transactions on Signal Processing], IEEE Transactions on , vol.23, no.1, pp.
67-72, Feb 1975. URL: \url{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1162641} 
 
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

print.stepPattern(symmetric2)



## Same example seen in dtw help, but with asymmetric  step pattern

idx<-seq(0,6.28,len=100);
query<-sin(idx)+runif(100)/10;
template<-cos(idx);

## Do the computation (takes ~5 s)
asy<-dtw(query,template,keep=TRUE,step=asymmetric);

dtwPlot(asy,type="density",main="Sine and cosine, asymmetric step")





}




\concept{Dynamic Time Warp}
\concept{Dynamic Programming}
\concept{Step pattern}
\concept{Local constraint}
\concept{Asymmetric DTW}
\concept{Symmetric DTW}



\keyword{ ts }
    
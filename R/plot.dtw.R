
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




#' Plotting of dynamic time warp results
#' 
#' Methods for plotting dynamic time warp alignment objects returned by
#' [dtw()].
#' 
#' `dtwPlot` displays alignment contained in `dtw` objects.
#' 
#' Various plotting styles are available, passing strings to the `type`
#' argument (may be abbreviated):
#' 
#'  * `alignment` plots the warping curve in `d`;
#'  * `twoway` plots a point-by-point comparison, with matching lines; see [dtwPlotTwoWay()];
#'  * `threeway` vis-a-vis inspection of the timeseries and their warping curve; see [dtwPlotThreeWay()];
#'  * `density` displays the cumulative cost landscape with the warping path overimposed; see [dtwPlotDensity()]
#' 
#' Additional parameters are passed to the plotting functions: use with
#' care.
#' 
#' @aliases dtwPlot dtwPlotAlignment  plot.dtw 
#' @param x,d `dtw` object, usually result of call to [dtw()]
#' @param xlab label for the query axis
#' @param ylab label for the reference axis
#' @param type general style for the plot, see below
#' @param plot.type type of line to be drawn, used as the `type` argument
#' in the underlying `plot` call
#' @param ... additional arguments, passed to plotting functions
#' @section Warning: These functions are incompatible with mechanisms for
#' arranging plots on a device: `par(mfrow)`, `layout` and
#' `split.screen`.
#' @author Toni Giorgino
#' @family plot
#' @seealso [dtwPlotTwoWay()], [dtwPlotThreeWay()] and  [dtwPlotDensity()] for details 
#' on the respective plotting styles.
#' @keywords ts hplot
#' @examples
#' 
#' ## Same example as in dtw
#' 
#' idx<-seq(0,6.28,len=100);
#' query<-sin(idx)+runif(100)/10;
#' reference<-cos(idx)
#' 
#' alignment<-dtw(query,reference,keep=TRUE);
#' 
#' # A sample of the plot styles. See individual plotting functions for details
#' 
#' plot(alignment, type="alignment",
#'   main="DTW sine/cosine: simple alignment plot")
#'   
#' plot(alignment, type="twoway",
#'   main="DTW sine/cosine: dtwPlotTwoWay")
#' 
#' plot(alignment, type="threeway",
#'   main="DTW sine/cosine: dtwPlotThreeWay")
#'   
#' plot(alignment, type="density",
#'   main="DTW sine/cosine: dtwPlotDensity")
#'   
#' @name dtwPlot
#' @import graphics
#' @export
plot.dtw <- function(x, type="alignment", ...) {
  
  pt<-pmatch(type,c("alignment",
                    "twoway",
                    "threeway",
                    "density"));
  switch(pt, 	dtwPlotAlignment(x, ...),
                dtwPlotTwoWay(x, ...),
                dtwPlotThreeWay(x, ...),
                dtwPlotDensity(x, ...)
 );
}


## an alias
#' @export
dtwPlot <- plot.dtw;


#' @rdname dtwPlot
#' @export 
dtwPlotAlignment <- function(d, xlab="Query index", ylab="Reference index", plot.type="l", ...) {
  plot( d$index1,d$index2,
        xlim=c(1,d$N),ylim=c(1,d$M),
	xlab=xlab,ylab=ylab,type=plot.type,
        ...
	);
}



## Well-known and much-copied pairwise matching


#' Plotting of dynamic time warp results: pointwise comparison
#' 
#' Display the query and reference time series and their alignment, arranged
#' for visual inspection.
#' 
#' The two vectors are displayed via the [matplot()] functions; their
#' appearance can be customized via the `type` and `pch` arguments
#' (constants or vectors of two elements).  If `offset` is set, the
#' reference is shifted vertically by the given amount; this will be reflected
#' by the *right-hand* axis.
#' 
#' Argument `match.indices` is used to draw a visual guide to matches; if
#' a vector is given, guides are drawn for the corresponding indices in the
#' warping curve (match lines). If integer, it is used as the number of guides
#' to be plotted. The corresponding style is customized via the
#' `match.col` and `match.lty` arguments.
#' 
#' If `xts` and `yts` are not supplied, they will be recovered from
#' `d`, as long as it was created with the two-argument call of
#' [dtw()] with `keep.internals=TRUE`.  Only single-variate time
#' series can be plotted this way.
#' 
#' @param d an alignment result, object of class `dtw`
#' @param xts query vector
#' @param yts reference vector
#' @param xlab,ylab axis labels
#' @param offset displacement between the timeseries, summed to reference
#' @param match.col,match.lty color and line type of the match guide lines
#' @param match.indices indices for which to draw a visual guide
#' @param ts.type,pch graphical parameters for timeseries plotting, passed to
#' `matplot`
#' @param ... additional arguments, passed to `matplot`
#' @note When `offset` is set values on the left axis only apply to the
#' query.
#' @section Warning: The function is incompatible with mechanisms for arranging
#' plots on a device: `par(mfrow)`, `layout` and `split.screen`.
#' @author Toni Giorgino
#' @family plot
#' @seealso [dtwPlot()] for other dtw plotting functions,
#' [matplot()] for graphical parameters.
#' @keywords hplot
#' @examples
#' 
#' 
#' ## A noisy sine wave as query
#' ## A cosine is for reference; sin and cos are offset by 25 samples
#' 
#' idx<-seq(0,6.28,len=100);
#' query<-sin(idx)+runif(100)/10;
#' reference<-cos(idx)
#' dtw(query,reference,step=asymmetricP1,keep=TRUE)->alignment;
#' 
#' 
#' ## Equivalent to plot(alignment,type="two");
#' dtwPlotTwoWay(alignment);
#' 
#' 
#' ## Highlight matches of chosen QUERY indices. We will do some index
#' ## arithmetics to recover the corresponding indices along the warping
#' ## curve
#' 
#' hq <- (0:8)/8              
#' hq <- round(hq*100)      #  indices in query for  pi/4 .. 7/4 pi
#' 
#' hw <- (alignment$index1 %in% hq)   # where are they on the w. curve?
#' hi <- (1:length(alignment$index1))[hw];   # get the indices of TRUE elems
#' 
#' 
#' ## Beware of the reference's y axis, may be confusing
#' plot(alignment,offset=-2,type="two", lwd=3, match.col="grey50",
#'      match.indices=hi,main="Match lines shown every pi/4 on query");
#' 
#' legend("topright",c("Query","Reference (rt. axis)"), pch=21, col=1:6)
#' 
#' 
#' 
#' @export 
dtwPlotTwoWay <- function(d,xts=NULL,yts=NULL, offset=0,
			ts.type="l",pch=21, 
                        match.indices=NULL,
			match.col="gray70", match.lty=3,
			xlab="Index", ylab="Query value", 
			... ) {

	if(is.null(xts) || is.null(yts))  {
            xts <- d$query;
            yts <- d$reference;
        }
    
	if(is.null(xts) || is.null(yts)) 
		stop("Original timeseries are required");

        ytso<-yts+offset;

        ## pad to longest
        maxlen<-max(length(xts),length(ytso));
        length(xts)<-maxlen;
        length(ytso)<-maxlen;

	
	## save default, for resetting...
	def.par <- par(no.readonly = TRUE);

        ## make room for secondary axis, if any
        if(offset!=0) {
          par(mar=c(5,4,4,4)+.1);
        }

	## plot q+t
	matplot(cbind(xts,ytso),
                type=ts.type,pch=pch, 
                xlab=xlab, ylab=ylab,
                axes=FALSE,
                ...);

        ## box and main axis
        ## compute range covering all values
        box();
        axis(1);
        axis(2,at=pretty(xts));

        ## display secondary axis if offset
        if(offset!=0) {
          rightTicks <- pretty(yts);
          axis(4,at=rightTicks+offset,labels=rightTicks);
        }


	## plot the matching 
	# par(par.match);
        if(is.null(match.indices)) {
          ml<-length(d$index1);
          idx<-1:ml;
        } else if(length(match.indices)==1) {
          idx <- seq(from=1,
                     to=length(d$index1),
                     length.out=match.indices);
        } else {
          idx <- match.indices;
        }

	## x0, y0 	coordinates of points from which to draw.
	## x1, y1 	coordinates of points to which to draw.
	segments(d$index1[idx],xts[d$index1[idx]],
		 d$index2[idx],ytso[d$index2[idx]],
		 col=match.col,lty=match.lty);

	
	par(def.par)#- reset to default

}







## ##################################################
## Global distance density plot

# for each plot, we should set: color, width, style, type
# for match lines: color, width, style


#' Plotting of dynamic time warp results: annotated warping function
#' 
#' Display the query and reference time series and their warping curve,
#' arranged for visual inspection.
#' 
#' 
#' The query time series is plotted in the bottom panel, with indices growing
#' rightwards and values upwards. Reference is in the left panel, indices
#' growing upwards and values leftwards. The warping curve panel matches
#' indices, and therefore element (1,1) will be at the lower left, (N,M) at the
#' upper right.
#' 
#' Argument `match.indices` is used to draw a visual guide to matches; if
#' a vector is given, guides are drawn for the corresponding indices in the
#' warping curve (match lines). If integer, it is used as the number of guides
#' to be plotted. The corresponding style is customized via the
#' `match.col` and `match.lty` arguments.
#' 
#' If `xts` and `yts` are not supplied, they will be recovered from
#' `d`, as long as it was created with the two-argument call of
#' [dtw()] with `keep.internals=TRUE`.  Only single-variate time
#' series can be plotted.
#' 
#' @param d an alignment result, object of class `dtw`
#' @param xts query vector
#' @param yts reference vector
#' @param xlab label for the query axis
#' @param ylab label for the reference axis
#' @param main main title
#' @param type.align line style for warping curve plot
#' @param type.ts line style for timeseries plot
#' @param match.indices indices for which to draw a visual guide
#' @param margin outer figure margin
#' @param inner.margin inner figure margin
#' @param title.margin space on the top of figure
#' @param ... additional arguments, used for the warping curve
#' @family plot
#' @section Warning: The function is incompatible with mechanisms for arranging
#' plots on a device: `par(mfrow)`, `layout` and `split.screen`.
#' Appearance of the match lines and timeseries currently can not be
#' customized.
#' @author Toni Giorgino
#' @keywords hplot
#' @examples
#' 
#' 
#' ## A noisy sine wave as query
#' ## A cosine is for reference; sin and cos are offset by 25 samples
#' 
#' idx<-seq(0,6.28,len=100);
#' query<-sin(idx)+runif(100)/10;
#' reference<-cos(idx)
#' dtw(query,reference,keep=TRUE)->alignment;
#' 
#' 
#' ## Beware of the reference's y axis, may be confusing
#' ## Equivalent to plot(alignment,type="three");
#' dtwPlotThreeWay(alignment);
#' 
#' 
#' ## Highlight matches of chosen QUERY indices. We will do some index
#' ## arithmetics to recover the corresponding indices along the warping
#' ## curve
#' 
#' hq <- (0:8)/8              
#' hq <- round(hq*100)      #  indices in query for  pi/4 .. 7/4 pi
#' 
#' hw <- (alignment$index1 %in% hq)   # where are they on the w. curve?
#' hi <- (1:length(alignment$index1))[hw];   # get the indices of TRUE elems
#' 
#' dtwPlotThreeWay(alignment,match.indices=hi);
#' 
#' @export 
dtwPlotThreeWay <- function(d,xts=NULL,yts=NULL,
                            type.align="l",type.ts="l",
                            match.indices=NULL,
                            margin=4, inner.margin=0.2, title.margin=1.5,
                            xlab="Query index",ylab="Reference index",main="Timeseries alignment",
                            ... ) {

     if(is.null(xts) || is.null(yts))  {
         xts <- d$query;
         yts <- d$reference;
     }

     # Sanity check
     if(is.null(xts) || is.null(yts))
       stop("Original timeseries are required");

     # Coerce to plain vectors
     xts <- as.matrix(xts);
     yts <- as.matrix(yts);

     # Verify if not multivariate
     if( ncol(xts)>1 || ncol(yts)>1 )
       stop("Only single-variate timeseries can be displayed. (You may want to extract a column for visualization purposes.)");


     def.par <- par(no.readonly = TRUE) # save default, for resetting...

     layout(matrix(c(3,1,0,2),2,2,byrow=TRUE), c(1,3), c(3,1), TRUE);

     
     imar<-inner.margin;
     
     bmar<-margin;
     lmar<-margin;
     tmar<-margin+title.margin;
     rmar<-margin;

     mlab=margin/2;
     mtex=margin/6;

     nn<-length(xts);
     mm<-length(yts);

     
     # Plot the warping function
     par(mar=c(imar,imar,tmar,rmar));

     # todo: plot over segments

     plot(d$index1,d$index2,type=type.align,
          xlim=c(1,nn),ylim=c(1,mm),
          ax=FALSE,main=main, ...
          ); # fake a diagonal, to set the axes


     # vertical match segments
     #  1 value: plot total of N elements
     if(length(match.indices)==1) {
       match.indices <- seq(from=1,
                            to=length(d$index1),
                            length.out=match.indices);
     }

     #  vector: use specified indices
     if(! is.null(match.indices) ) {      # vertical match segments
       idx <- match.indices;
       segments(d$index1[idx],0,
                d$index1[idx],d$index2[idx],
                col="grey60",lty=3);
                                        # horz.
       segments(0,d$index2[idx],
                d$index1[idx],d$index2[idx],
                col="grey60",lty=3);
     }

     
     box();

     
     # axis are 1- bot; 2- left; 3- top; 4- right
     # Plot query (horizontal, bottom)
     par(mar=c(bmar,imar,imar,rmar));

     plot(xts ~ c(1:nn), type=type.ts,
          xlab=xlab ,mgp=c(mlab,mtex,0) ,ax=FALSE,
          );
     axis(1);
     axis(2);
     box();

     # Plot reference (vertical, left)
     par(mar=c(imar,lmar,tmar,imar));

     # reverse the horiz. axis so that rotation is more natural
     plot(c(1:mm) ~ yts, xlim=rev(range(yts)), type=type.ts,
          ylab=ylab, mgp=c(mlab,mtex,0) , ax=FALSE,
          );
     axis(3);
     axis(2);
     box();

     par(def.par)#- reset to default

}







#' Display the cumulative cost density with the warping path overimposed
#' 
#' The plot is based on the cumulative cost matrix. It displays the optimal alignment
#' as a "ridge" in the global cost landscape.
#' 
#' The alignment must have been
#' constructed with the `keep.internals=TRUE` parameter set. 
#' 
#' If `normalize` is `TRUE`, the *average* cost per step is
#' plotted instead of the cumulative one. Step averaging depends on the
#' [stepPattern()] used.
#' 
#' @param d an alignment result, object of class `dtw`
#' @param normalize show per-step average cost instead of cumulative cost
#' @param xlab label for the query axis
#' @param ylab label for the reference axis
#' @param ... additional parameters forwarded to plotting functions
#' @family plot
#' @examples 
#' 
#' ## A study of the "Itakura" parallelogram
#' ##
#' ## A widely held misconception is that the "Itakura parallelogram" (as
#' ## described in the original article) is a global constraint.  Instead,
#' ## it arises from local slope restrictions. Anyway, an "itakuraWindow",
#' ## is provided in this package. A comparison between the two follows.
#' 
#' ## The local constraint: three sides of the parallelogram are seen
#' idx<-seq(0,6.28,len=100);
#' query<-sin(idx)+runif(100)/10;
#' reference<-cos(idx)
#' 
#' ita <- dtw(query,reference,keep=TRUE,step=typeIIIc)
#' dtwPlotDensity(ita, main="Slope-limited asymmetric step (Itakura)")
#' 
#' ## Symmetric step with global parallelogram-shaped constraint. Note how
#' ## long (>2 steps) horizontal stretches are allowed within the window.
#' 
#' dtw(query,reference,keep=TRUE,window=itakuraWindow)->ita;
#' dtwPlotDensity(ita,
#'         main="Symmetric step with Itakura parallelogram window")
#' 
#' @export 
dtwPlotDensity <- function(d, normalize=FALSE,
                           xlab="Query index", ylab="Reference index", ...) {
    
    cm<-d$costMatrix;
    
    if(is.null(cm)) 
        stop("dtwPlotDensity requires dtw internals (set keep.internals=TRUE on dtw() call)");
    
    ## We can safely modify cm locally
    if(normalize) {
        norm <- attr(d$stepPattern,"norm");
        if(is.na(norm))
            stop("No normalization known for step pattern used");
        
        if(norm=="N") {
            cm <- cm / row(cm);
        } else if(norm=="N+M") {
            cm <- cm / (row(cm)+col(cm));
        } else if(norm=="M") {
            cm <- cm / col(cm);
        }
    }
    
    xd<-dim(cm)[1];
    yd<-dim(cm)[2];
    
    image(cm,col=grDevices::terrain.colors(100),x=1:xd,y=1:yd,
          xlab=xlab,ylab=ylab, ...);
    contour(cm,x=1:xd,y=1:yd,add=TRUE);
    lines(d$index1,d$index2,col="blue",lwd=2);
}




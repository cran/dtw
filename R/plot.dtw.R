###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino@gmail.com>           #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: plot.dtw.R 128 2008-05-29 14:00:10Z tonig $
#                                                             #
###############################################################




## Plot a dtw non-object, switching depending on requested type

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
dtwPlot <- plot.dtw;




dtwPlotAlignment <- function(d, xlab="Query index", ylab="Template index", ...) {
  plot( d$index1,d$index2,
        xlim=c(1,d$N),ylim=c(1,d$M),
	xlab=xlab,ylab=ylab, ...
	);
}


## Normalization plots the average cost per step instead of
## the cumulative cost

dtwPlotDensity <- function(d, normalize="no",
                           xlab="Query index", ylab="Template index", ...) {

    cm<-d$costMatrix;

    if(is.null(cm)) 
      stop("dtwPlotDensity requires dtw internals (set keep.internals=TRUE on dtw() call)");

    ## We can safely modify cm locally
    nt<-pmatch(normalize,c("N","N+M"));

    if(!is.na(nt)) {
        cm <- cm / switch(nt,
                          row(cm),
                          row(cm)+col(cm));
    }

    xd<-dim(cm)[1];
    yd<-dim(cm)[2];

    image(cm,col=terrain.colors(100),x=1:xd,y=1:yd,
          xlab=xlab,ylab=ylab, ...);
    contour(cm,x=1:xd,y=1:yd,add=TRUE);
    lines(d$index1,d$index2,col="blue",lwd=2);
}




## Well-known and much-copied pairwise matching

dtwPlotTwoWay <- function(d,xts=NULL,yts=NULL, offset=0,
			type="o",pch=21, 
			xlab="Index", ylab="Query value", 
			match.col="gray70",
			... ) {

	if(is.null(xts) || is.null(yts))
		stop("Original timeseries are required");
        ytso<-yts+offset;

        ## pad to longest
        maxlen<-max(length(xts),length(ytso));
        length(xts)<-maxlen;
        length(ytso)<-maxlen;
	
	## save default, for resetting...
	def.par <- par(no.readonly = TRUE);


	## plot q+t
	matplot(cbind(xts,ytso),
			type=type,pch=pch, 
			xlab=xlab, ylab=ylab,
			...);

        ## display secondary axis if offset
        if(offset!=0) {
          axt<-axTicks(2);
          axis(4,at=axt+offset,labels=axt);
        }


	## plot the matching 
	# par(par.match);
	ml<-length(d$index1);
	idx<-1:ml;

	## x0, y0 	coordinates of points from which to draw.
	## x1, y1 	coordinates of points to which to draw.
	segments(d$index1[idx],xts[d$index1[idx]],
		 d$index2[idx],ytso[d$index2[idx]],
		 col=match.col);

	
	par(def.par)#- reset to default

}







## ##################################################
## Global distance density plot


dtwPlotThreeWay <- function(d,xts=NULL,yts=NULL,type.align="p",type.ts="l",
			margin=4, inner.margin=0.2, title.margin=1.5,
			xlab="Query index",ylab="Template index",main="Timeseries alignment",
			... ) {

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
     plot(d$index1,d$index2,,type=type.align,
          xlim=c(1,nn),ylim=c(1,mm),
          ax=FALSE,main=main); # fake a diagonal, to set the axes
     # lines(d$index1,d$index2);
     box();

     # axis are 1- bot; 2- left; 3- top; 4- right
     # Plot query (horizontal, bottom)
     par(mar=c(bmar,imar,imar,rmar));
     plot(xts ~ c(1:nn) , type=type.ts ,xlab=xlab ,mgp=c(mlab,mtex,0) ,ax=FALSE );
     axis(1);
     axis(2);
     box();

     # Plot template (vertical, left)
     par(mar=c(imar,lmar,tmar,imar));
     # reverse the horiz. axis so that rotation is more natural
     plot(c(1:mm) ~ yts , xlim=rev(range(yts)), type=type.ts,
          ylab=ylab, mgp=c(mlab,mtex,0) , ax=FALSE );
     axis(3);
     axis(2);
     box();

     par(def.par)#- reset to default

}


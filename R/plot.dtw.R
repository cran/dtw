###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino@gmail.com>           #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: plot.dtw.R 51 2007-12-11 10:59:08Z tonig $
#                                                             #
###############################################################




## Plot a dtw non-object, switching depending on requested type

plot.dtw <- function(x, type="alignment", ...) {
  
  pt<-pmatch(type,c("alignment",
                    "twoway",
                    "threeway",
                    "density"));
  switch(pt, 	dtwPlotAlignment(x, ...),
                .dtwPlotTwoWay(x, ...),
                dtwPlotThreeWay(x, ...),
                dtwPlotDensity(x, ...)
 );
}


## an alias
dtwPlot <- plot.dtw;




dtwPlotAlignment <- function(d, xlab="Query index", ylab="Template index", ...) {
  plot( d$index1,d$index2,
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

	## We can safely modify x locally
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
	lines(d$index1,d$index2,col="red",lwd=2);
}


.dtwPlotTwoWay <- function(x) {
  stop("Not implemented");
}







## ##################################################
## Global distance density plot


dtwPlotThreeWay <- function(d,xts=NULL,yts=NULL,type.align="p",type.ts="l",
			margin=4, inner.margin=0.2, title.margin=1.5,
			xlab="Query index",ylab="Template index",main="Timeseries alignment",
			... ) {

#     xts<-x; yts<-y;
     if(is.null(xts) || is.null(yts))
       stop("Original timeseries are required");

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

     par(mar=c(imar,imar,tmar,rmar));
     plot(d$index1,d$index2,type=type.align,main=main,ax=FALSE);
     box();

     par(mar=c(bmar,imar,imar,rmar));
     plot(xts ~ c(1:nn) , type=type.ts ,xlab=xlab ,mgp=c(mlab,mtex,0) ,ax=FALSE );
     axis(1);
     axis(2);
     box();

     par(mar=c(imar,lmar,tmar,imar));
     plot(c(1:mm) ~ yts , type=type.ts ,ylab=ylab ,mgp=c(mlab,mtex,0) , ax=FALSE );
     axis(3);
     axis(2);
     box();

     par(def.par)#- reset to default

}


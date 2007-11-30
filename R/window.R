###############################################################
#                                                             #
#   Author: Toni Giorgino <toni.giorgino@gmail.com>           #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: window.R 51 2007-12-11 10:59:08Z tonig $
#                                                             #
###############################################################

##
## Document in dtwWindowingFunctions.Rd
##




## no warping window: no restrictions

`noWindow` <-
function(iw,jw,...) {
  return(TRUE);
}



## A band around the diagonal. The band includes the diagonal +-
## window.size.

`sakoeChibaWindow` <-
function(iw,jw,window.size,...) {
  return(abs(jw-iw)<=window.size);
}


## A band around the diagonal. The band includes the segment
## connecting [1,1] to [query.size,template.size] window.size,
## measured along the second index (columns)

`slantedBandWindow` <-
function(iw,jw,query.size,template.size,window.size,...) {
  diagj<-(iw*template.size/query.size);
  return(abs(jw-diagj)<=window.size);
}



## "Itakura" parallelogram: see documentation.
## 

`itakuraWindow` <- 
function(iw,jw,query.size,template.size,...) {
	## Shorthands
  	n<-query.size; 	
  	m<-template.size;

	ok<- 	(jw <  2*iw) &
 		(iw <= 2*jw) &
		(iw >= n-1-2*(m-jw)) &
		(jw >  m-1-2*(n-iw)) ;

	return(ok);
}


## Plot a sample of the windowing function

dtwWindow.plot <- function(fun,query.size=200,template.size=220,...) {
	n<-query.size;
  	m<-template.size;

	mm<-matrix(0,n,m);
	mm[fun(row(mm),col(mm),query.size=n,template.size=m,...)]<-1;

	image(  z=mm,
		x=1:n,
		y=1:m,
		xlab=sprintf("Query: samples 1..%d",n),
		ylab=sprintf("Template: samples 1..%d",m), 
		...
	 );
}



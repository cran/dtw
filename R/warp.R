###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino@gmail.com>           #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: plot.dtw.R 55 2007-12-13 12:25:49Z tonig $
#                                                             #
###############################################################


## Apply the warping curve to a given timeseries. If the curve is
## multi-valued, average multiple points, with a warning. Argument
## "inverse=T" warps the inverse of the curve (ie a template into a
## query).

## Direct:
##  all points in query have at least one image (NO?)
##  only one image, if asymmetric
##  if partial, image <= template
##  x should be as long as the range in ix
##  there could be gaps in ix, to be interpolated,
##  depending on the step pattern

## for each point in the template space, do an interpolated lookup
## into the y->x mapping

## Inverse: considering template as domain
##  if partial, some trailing points may have no image
##  otherwise, all points have one or more images


## sortedXyData

warp <- function(d,index.template=FALSE) {

  if(!is.dtw(d))
    stop("dtw object required");
  

  if(!index.template) {
    ## warp QUERY into template space
    iset<-d$index1;
    jset<-d$index2;
  } else {
    ## warp TEMPLATE into query
    iset<-d$index2;
    jset<-d$index1;
  }
  jmax<-max(jset);

  ## rebuild  index, interpolating holes
  ii<-approx(x=jset,y=iset,1:jmax);
  return(ii$y);    
     
}


## > ss<-warp.dtw(al,query)
## Warning message:
## In approx(x = jset, y = iset, 1:jmax) : collapsing to unique 'x' values
## > plot(template);lines(query)
## > plot(template);lines(ss)
## > st<-warp.dtw(al,template,inverse=T)
## Warning message:
## In approx(x = jset, y = iset, 1:jmax) : collapsing to unique 'x' values
## > Cairo()
## > plot(query);lines(st)
## > plot(query,type="l");points(st)
## > 

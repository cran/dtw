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


## Compute the (approximate) area between the warping function and the
## diagonal (in unit steps). 


warpArea <- function(d) {

  if(!is.dtw(d))
    stop("dtw object required");
  
  ## rebuild query->templ map, interpolating holes
  ii<-approx(x=d$index1,y=d$index2,1:d$N);

  dg<-seq(from=1,to=d$M,len=d$N);

  ad<-abs(ii$y-dg);
  sum(ad);
     
}

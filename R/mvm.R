###############################################################
#                                                             #
#   Author: Toni Giorgino <toni.giorgino,gmail.com>           #
#       Istituto di Neuroscienze (IN-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id: mvm.R 425 2016-08-25 19:48:58Z tonig $
#                                                             #
###############################################################



#############################
## Step pattern for MVM (minimum variance matching)
##

mvmStepPattern <- function(elasticity=20) {
  size <- elasticity;
  
  pn <- rep(1:size,each=2); # 112233...
  dx <- rep(c(1,0),size);   # 101010
  dy <- dx * ( (1:(2*size)) +1)/2;  # 102030
  w <-  rep(c(-1,1),size);
  
  tmp <- cbind(pn,dx,dy,w);
#  tmp <- as.vector(t(tmp));
  sp <- stepPattern(tmp);
  
  attr(sp,"norm") <- "N";
  attr(sp,"call") <- match.call();
  return(sp);
}




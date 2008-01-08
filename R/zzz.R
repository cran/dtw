###############################################################
#                                                             #
#   Author: Toni Giorgino <toni.giorgino@gmail.com>           #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: zzz.R 94 2008-01-08 16:44:46Z tonig $
#                                                             #
###############################################################

.First.lib <- function(lib, pkg)  {
  cat(paste("This is dtw ",
            utils::packageDescription("dtw")$Version,
	    ". Please see citation(\"dtw\") for use in publications.",
	    "\n",
            sep=""))
  library.dynam("dtw");
  library(proxy);

  ## Register DTW as a distance function into package proxy
  pr_DB$set_entry(FUN=dtwpairdist, names="DTW", description="Dynamic Time Warping",
                  loop=TRUE, formula="minimum of sum(x[xw[i]]-y[yw[i]]) over all xw, yw");

  invisible()
}

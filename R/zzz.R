###############################################################
#                                                             #
#   Author: Toni Giorgino <toni.giorgino@gmail.com>           #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: zzz.R 98 2008-01-09 14:20:18Z tonig $
#                                                             #
###############################################################

.First.lib <- function(lib, pkg)  {
  library(proxy);

  cat(paste("Loaded dtw v",
            utils::packageDescription("dtw")$Version,
	    ". See ?dtw for help, citation(\"dtw\") for use in publication.",
	    "\n",
            sep=""))
  library.dynam("dtw");

  ## Register DTW as a distance function into package proxy
  pr_DB$set_entry(FUN=dtwpairdist, names="DTW", 
                  loop=TRUE, type="metric",
                  description="Dynamic Time Warping",
                  formula="minimum of sum(x[xw[i]]-y[yw[i]]) over all monotonic xw, yw");

  invisible()
}

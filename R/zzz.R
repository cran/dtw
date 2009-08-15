###############################################################
#                                                             #
#   Author: Toni Giorgino <toni.giorgino@gmail.com>           #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: zzz.R 235 2009-06-30 15:13:52Z tonig $
#                                                             #
###############################################################

.onLoad <- function(lib, pkg)  {
  # library(proxy);

  cat(sprintf("Loaded dtw v%s. See ?dtw for help, citation(\"dtw\") for usage conditions.\n",
            utils::packageDescription("dtw")$Version ) );
      
  ## Register DTW as a distance function into package proxy
  pr_DB$set_entry(FUN=dtwpairdist, names="DTW", 
                  loop=TRUE, type="metric",
                  description="Dynamic Time Warping",
                  formula="minimum of sum(x[xw[i]]-y[yw[i]]) over all monotonic xw, yw");

  invisible()
}

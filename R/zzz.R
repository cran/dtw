###############################################################
#                                                             #
#   Author: Toni Giorgino <toni.giorgino,gmail.com>           #
#       Istituto di Ingegneria Biomedica (ISIB-CNR)           #
#       Consiglio Nazionale delle Ricerche                   #
#       www.isib.cnr.it                                       #
#                                                             #
#   $Id: zzz.R 267 2012-08-12 14:37:26Z tonig $
#                                                             #
###############################################################

.onAttach <- function(lib, pkg)  {

  packageStartupMessage(sprintf("Loaded dtw v%s. See ?dtw for help, citation(\"dtw\") for use in publication.\n",
            utils::packageDescription("dtw")$Version ) );
      
  ## Register DTW as a distance function into package proxy
  pr_DB$set_entry(FUN=dtwpairdist, names="DTW", 
                  loop=TRUE, type="metric",
                  description="Dynamic Time Warping",
                  formula="minimum of sum(x[xw[i]]-y[yw[i]]) over all monotonic xw, yw");

  # invisible()
}

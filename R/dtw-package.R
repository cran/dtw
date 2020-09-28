
##
## Copyright (c) 2006-2019 of Toni Giorgino
##
## This file is part of the DTW package.
##
## DTW is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## DTW is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
## or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
## License for more details.
##
## You should have received a copy of the GNU General Public License
## along with DTW.  If not, see <http://www.gnu.org/licenses/>.
##


#' Comprehensive implementation of Dynamic Time Warping (DTW) algorithms in R.
#' 
#' The DTW algorithm computes the stretch of the time axis which optimally maps
#' one given timeseries (query) onto whole or part of another (reference). It
#' yields the remaining cumulative distance after the alignment and the
#' point-by-point correspondence (warping function). DTW is widely used e.g.
#' for classification and clustering tasks in econometrics, chemometrics and
#' general timeseries mining.
#' 
#' Please see documentation for function [dtw()], which is the main
#' entry point to the package.
#' 
#' The R implementation in dtw provides:
#' 
#'  * arbitrary windowing functions (global constraints), eg. the
#' Sakoe-Chiba band; see [dtwWindowingFunctions()] 
#'  * arbitrary
#' transition types (also known as step patterns, slope constraints, local
#' constraints, or DP-recursion rules). This includes dozens of well-known
#' types; see [stepPattern()]: 
#'    * all step patterns classified by Rabiner-Juang, Sakoe-Chiba, and Rabiner-Myers; 
#'    * symmetric and asymmetric; 
#'    * Rabiner's smoothed variants; 
#'    * arbitrary, user-defined slope constraints 
#'  * partial matches: open-begin, open-end, substring matches 
#'  * proper, pattern-dependent, normalization (exact average distance per step) 
#'  * the Minimum Variance Matching (MVM) algorithm (Latecki et al.) 
#' 
#' Multivariate timeseries can be aligned with arbitrary local distance
#' definitions, leveraging the [proxy::dist()] function of package
#' \CRANpkg{proxy}. DTW itself becomes a distance function with the dist semantics.
#' 
#' In addition to computing alignments, the package provides:
#'  * methods for plotting alignments and warping functions in several classic
#' styles (see plot gallery); 
#'  * graphical representation of step patterns;
#'  * functions for applying a warping function, either direct or inverse;
#' and more. 
#' 
#' If you use this software, please cite it according to
#' `citation("dtw")`.  The package home page is at
#' <https://dynamictimewarping.github.io>.
#' 
#' @name dtw-package
#' @docType package
#' @author Toni Giorgino
#'
#' 
#' @seealso [dtw()] for the main entry point to the package;
#' [dtwWindowingFunctions()] for global constraints;
#' [stepPattern()] for local constraints;
#' [proxy::dist()], `analogue::distance()`, `vegan::vegdist()` to build local
#' cost matrices for multivariate timeseries and custom distance functions.
#' @references
#'  * Toni Giorgino. *Computing and Visualizing Dynamic Time
#' Warping Alignments in R: The dtw Package.* Journal of Statistical Software,
#' 31(7), 1-24. \doi{10.18637/jss.v031.i07}
#'  * Tormene, P.; Giorgino, T.; Quaglini, S. & Stefanelli, M. *Matching incomplete time
#' series with dynamic time warping: an algorithm and an application to
#' post-stroke rehabilitation.* Artif Intell Med, 2009, 45, 11-34 
#' \doi{10.1016/j.artmed.2008.11.007}
#'  * Rabiner, L. R., & Juang, B.-H. (1993). Chapter 4 in *Fundamentals of
#' speech recognition.* Englewood Cliffs, NJ: Prentice Hall.
#' @keywords package ts
#' @import proxy
#' @examples
#' 
#'  library(dtw);
#'  ## demo(dtw);
#' 
NULL



.onAttach <- function(lib, pkg)  {
    
    packageStartupMessage(sprintf("Loaded dtw v%s. See ?dtw for help, citation(\"dtw\") for use in publication.\n",
                                  utils::packageDescription("dtw")$Version ) );
    
    ## Register DTW as a distance function into package proxy
    if(!proxy::pr_DB$entry_exists("DTW")) {
        proxy::pr_DB$set_entry(FUN=dtwpairdist, names=c("DTW","dtw"), 
                               loop=TRUE, type="metric",
                               description="Dynamic Time Warping",
                               reference="Giorgino T (2009). Computing and Visualizing Dynamic Time Warping Alignments in R: The dtw Package. Journal of Statistical Software, 31 (7), pp. 1--24. <URL: http://www.jstatsoft.org/v31/i07/>.",
                               formula="minimum of sum(x[xw[i]]-y[yw[i]]) over all monotonic xw, yw");
    }
    
    # invisible()
}






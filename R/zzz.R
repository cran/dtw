.First.lib <- function(lib, pkg)  {
  cat(paste("This is dtw ",
            utils::packageDescription("dtw")$Version,
	    ". Please see citation(\"dtw\") on how to cite in publications.",
	    "\n",
            sep=""))
  library.dynam("dtw");
  invisible()
}

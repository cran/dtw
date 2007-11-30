.First.lib <- function(lib, pkg)  {
  cat(paste("This is dtw",
            utils::packageDescription("dtw")$Version,
            "\n"))
  library.dynam("dtw");
  invisible()
}

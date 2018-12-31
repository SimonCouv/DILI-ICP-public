get_versioned_fig_dir <- function(figname, ext, filepath=getwd()){
  library(readr)
  fullname <- paste0(figname, format(Sys.time(), "_v%Y%m%d_%H%M%S"), ".",ext)
  file.path(filepath,fullname)
}
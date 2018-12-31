get_versioned_png_dir <- function(figname, filepath=getwd()){
  library(readr)
  fullname <- paste0(figname, format(Sys.time(), "_v%Y%m%d_%H%M%S"), ".png")
  file.path(filepath,fullname)
}
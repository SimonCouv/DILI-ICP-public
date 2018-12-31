save_versioned_data <- function(object, folder=getwd()){
  objectname <- deparse(substitute(object))
  fullname <- paste0(objectname, format(Sys.time(), "_v%Y%m%d_%H%M%S"), ".RDS")
  filedir <- file.path(folder,fullname)
  saveRDS(object, filedir)
  print(sprintf("This is save_versioned_data. Saved succesfully to %s", filedir))
}
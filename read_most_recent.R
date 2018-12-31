read_most_recent <- function(pattern, readdir=getwd()){
  library(stringr)
  #depends on versioning structure in save_versioned_data
  full.file.list <- list.files(path = readdir, pattern = pattern, full.names = TRUE)
  matched <- str_match(full.file.list, "([^/]*)_v(\\d{8}_\\d{6}).*")
  if (nrow(matched)==0)
    stop(sprintf("No matches for pattern %s in directory %s.\n Directory %s contains:\n %s", pattern, readdir, readdir, paste0(dir(readdir),collapse="\n")))
  full.file.df <- data.frame(core_name = matched[,2], date_time = matched[,3])
  full.file.df$full_name <- full.file.list
  full.file.df$date_time <- as.POSIXct(full.file.df$date_time, format = "%Y%m%d_%H%M%S")
  to.read <- full.file.df %>% dplyr::group_by(core_name) %>% .[order(.$date_time, decreasing = TRUE),] %>% dplyr::slice(1)
  for (i in 1:nrow(to.read)){
    assign(as.character(to.read$core_name[i]),readRDS(file = to.read$full_name[i]), .GlobalEnv)
    print(sprintf("This is read_most_recent. Successfully read %s to object %s", to.read$full_name[i], to.read$core_name[i]))
  }
}
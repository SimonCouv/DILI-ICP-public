classify_chromHMM15_enh <- function(chrom15){
  #Ernst2015 fig6c
  if (any(!(chrom15 %in% c(0,4:7,12)))) 
    stop(sprintf("Classification not specified for state(s)"), paste0(unique(chrom15[!(chrom15 %in% c(0,4:7,12))]),collapse=","))
  ifelse(chrom15 %in% 4:5,
         "enh_strong",
         ifelse(chrom15 %in% 6:7,
                "enh_weak",
                NA))
}
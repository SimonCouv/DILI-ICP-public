classify_chromHMM15_prom <- function(chrom15){
  #Ernst2015 fig6c
  if (any(!(chrom15 %in% c(0,1:3,10)))) 
    stop(sprintf("Classification not specified for state(s)"), paste0(unique(chrom15[!(chrom15 %in% c(0,1:3,10))]),collapse=","))
  ifelse(chrom15 %in% 1:2,
         "prom",
         ifelse(chrom15 ==3,
                "poised",
                NA))
}
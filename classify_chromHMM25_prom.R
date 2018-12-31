classify_chromHMM25_prom <- function(chrom25){
  #Ernst2015 fig6c
  if (any(!(chrom25 %in% c(0:4,9,22,23)))) 
    stop(sprintf("Classification not specified for state(s)"), paste0(unique(chrom25[!(chrom25 %in% c(0:5,9,22,23))]),collapse=","))
  ifelse(chrom25==1,
         "TSS",
         ifelse(chrom25 %in% 2:4,
                "prom",
                ifelse(chrom25 %in% 22:23,
                       "poised_bivalent",
                       ifelse(chrom25==9,
                              "prom_enh_tx",
                              NA))))
}
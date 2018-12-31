classify_chromHMM25_enh <- function(chrom25){
  #Ernst2015 fig6c
  if (any(!(chrom25 %in% c(0,9:18)))) 
    stop(sprintf("Classification not specified for state(s)"), paste0(unique(chrom25[!(chrom25 %in% c(0,9:18))]),collapse=","))
  ifelse(chrom25==9,
         "prom_enh_tx",
         ifelse(chrom25 %in% 10:12,
                "enh_tx",
                ifelse(chrom25 %in% 13:15, 
                       "enh_strong", 
                       ifelse(chrom25 %in% 16:18, 
                              "enh_weak",
                              NA))))
}
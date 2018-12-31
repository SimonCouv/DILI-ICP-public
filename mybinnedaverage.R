mybinnedaverage <- function(bins, signal, varname="binned_score"){   #based on line 989 onw of https://github.com/jotsetung/Gviz/blob/master/R/Gviz-methods.R
  # if (!identical(seqlevels(bins), seqnames(signal)))
  #   stop("'seqlevels(bin)' and 'seqnames(numvar)' must be identical")
  sc <- values(signal)
  ol <- as.matrix(findOverlaps(bins, signal))
  scn <- sapply(split(ol[,2], ol[,1]), function(i) mean(sc[i,]), USE.NAMES=FALSE)
  sc <- matrix(NA, nrow = length(bins))
  sc[unique(ol[,1])] <- scn
  scdf <- DataFrame(binscore = sc)
  names(scdf) <- varname
  values(bins) <- scdf
  bins
}
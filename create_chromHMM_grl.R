create_chromHMM_grl <- function(gene.gr, prom.gr, enh.gr){
  gene.grl <- split(gene.gr, mcols(gene.gr))  #mcol contains gene symbol
  # gr<- append(subsetByOverlaps(chromHMM_25_prom.gr, gene.gr),
  #             subsetByOverlaps(chromHMM_25_enh.gr, gene.gr))
  grl <- list()
  for (genename in names(gene.grl)){
    gene.gr <- gene.grl[[genename]]
    gr <- subsetByOverlaps(append(prom.gr,enh.gr), gene.gr)
    if (length(reduce(gr)) > length(gr))
      stop(sprintf("found overlap in prom/enh ranges for gene %s", genename))
    else
      grl[[genename]] <- gr
  }
  grl
}
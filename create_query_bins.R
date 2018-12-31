create_query_bins <- function(query.gr, metacol ,flt, tilew){
  #query ranges, = Dixon ICP genes with >=1 DILI hit
  
  filtered.query.gr <- query.gr[elementMetadata(query.gr)[,metacol] %in% flt]
  seqlevels(filtered.query.gr) <- seqlevelsInUse(filtered.query.gr)
  filtered.query.gr <- sort(filtered.query.gr)
  #https://github.com/hochwagenlab/hwglabr2/wiki/Apply-function-to-GRanges-scores-in-genomic-tiles
  #create bins
  # bins <- GenomicRanges::tileGenome(seqlengths(filtered.query.gr), tilewidth = tilew, cut.last.tile.in.chrom = TRUE)
  bins_local <- GenomicRanges::tile(filtered.query.gr, width = tilew)
  names(bins_local) <- elementMetadata(filtered.query.gr)[,metacol]
  return(bins_local)
}
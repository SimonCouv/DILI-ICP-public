annotate_search_range <- function(gene.start, gene.end, strand, prom.len){
  if (strand == 1){
    range.start = gene.start - prom.len
    range.end = gene.end
  } else if (strand == -1){
    range.start = gene.start
    range.end = gene.end + prom.len
  }
  list(range.start = range.start, range.end = range.end)
}
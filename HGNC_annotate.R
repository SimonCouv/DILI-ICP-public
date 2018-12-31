HGNC_annotate <- function(genelist){
  
  for (i in 1:nrow(genelist)){
    print(i)
    # browser()
    dixongene <- genelist$gene[i]
    if (any(HGNC.table$symbol == dixongene) && sum(HGNC.table$symbol == dixongene) == 1){
      HGNC <- HGNC.table$hgnc_id[HGNC.table$symbol==dixongene]
      genelist[i, 'HGNC'] <- HGNC
    } else if (any(HGNC.table$symbol == dixongene) && sum(HGNC.table$symbol == dixongene)  > 1){
      stop(paste0("Multiple HGNC hits for input", dixongene))
    } else if (any(grepl(pattern = dixongene, x=HGNC.table$alias_symbol))){
      HGNC <- HGNC.table$alias_symbol[HGNC.table$symbol==dixongene]
      genelist[i, 'HGNC'] <- HGNC
    } else genelist[i, 'HGNC'] <- 'manual search'
  }
  genelist
}

in1KG <- Vectorize(function(rs){
  res <- entrez_search(db='snp', term = paste0(rs,'[RS] AND "1000genomes"[Submitter Handle]'))
  length(res$ids)>0
})
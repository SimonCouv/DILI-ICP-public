in1KG <- Vectorize(function(rs){
  res <- entrez_search(db='snp', term = paste0(rs,'[RS] AND "1000genomes"[Submitter Handle]'))
  length(res$ids)>0
})
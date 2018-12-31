lookup_rs <- function(rs.list,linktable){
  lookup.df <- expand.grid(
    SNP = rs.list,
    assembly = c("GRCh37", "GRCh38"),
    dataset = c("all", "MC", "HC", "AC")
  )
  
  match.list <- apply(lookup.df,1,search_GWAS_loci, linktable)
  # browser()
  lookup.df$match.rank <- sapply(match.list,function(x)x[[1]])
  lookup.df$match.offset <- sapply(match.list,function(x)x[[2]])
  lookup.df
}
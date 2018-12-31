get_BM_anno <- function(rs.list, ens37, ens38){
  
  BM.GRCh38 <- getBM(
    attributes = c("refsnp_id", "refsnp_source", "chr_name", "chrom_start", "chrom_end", "chrom_strand", "allele", "minor_allele", "minor_allele_freq", "associated_gene"),
    filters = c("snp_filter","chr_name"),
    values = list(rs.list, c(1:22, "X", "Y", "MT")),
    mart = ens38
  )
  
  BM.GRCh37 <- getBM(
    attributes = c("refsnp_id", "refsnp_source", "chr_name", "chrom_start", "chrom_strand", "allele", "minor_allele", "minor_allele_freq", "associated_gene"),
    filters = c("snp_filter","chr_name"),
    values = list(rs.list, c(1:22, "X", "Y", "MT")),
    mart = ens37
  )
  
  BM.GRCh37.s <- simplify_associated_genes(BM.GRCh37)
  BM.GRCh38.s <- simplify_associated_genes(BM.GRCh38)
  
  names(BM.GRCh37.s) <- paste0(names(BM.GRCh37.s),".GRCh37")
  names(BM.GRCh38.s) <- paste0(names(BM.GRCh38.s),".GRCh38")
  
  
  merge(
    BM.GRCh38.s, by.x = "refsnp_id.GRCh38",
    BM.GRCh37.s, by.y = "refsnp_id.GRCh37")
}
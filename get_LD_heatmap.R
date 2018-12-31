get_LD_heatmap <- function(genename, dili_icp_rs_ld_plink, SNP.gene.matches){
  
  gene_DILI_ld_long <- dili_icp_rs_ld_plink %>% 
    dplyr::filter(gene==genename & population =="CEU") %>% 
    dplyr::select(SNP_A, SNP_B, R2) %>%
    merge(.,SNP.gene.matches, by.x="SNP_A", by.y="RefSNP") %>%
    dplyr::filter(SNP_B %in% SNP.gene.matches$RefSNP[SNP.gene.matches$gene==genename])%>%
    dplyr::select(SNP_A, SNP_B, R2,chr_A=chr, bp_A=bp)
  
  rs_order <- order(unique(gene_DILI_ld_long$bp_A))
  rs_label_order <- unique(gene_DILI_ld_long$SNP_A)[rs_order]
  dist <- unique(gene_DILI_ld_long$bp_A)[rs_order]
  
  gene_DILI_ld_sq <- gene_DILI_ld_long %>%
    dplyr::select(SNP_A,SNP_B,R2)%>%
    distinct() %>%
    dcast(SNP_A~SNP_B, drop=FALSE, value.var="R2") %>%
    .[rs_order,rs_label_order] %>%
    as.matrix()
  rownames(gene_DILI_ld_sq) <- colnames(gene_DILI_ld_sq)
  
  return(LDheatmap(gene_DILI_ld_sq, flip = TRUE, genetic.distances = dist))
}
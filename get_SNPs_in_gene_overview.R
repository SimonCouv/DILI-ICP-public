get_SNPs_in_gene_overview <- function(genelist, ensSNP){
  
  find_SNP_in_gene_range <- function(row){
    gene <- row[1]
    assembly <- row[2]
    dataset <- row[3]
    
    
    SNP.table <- switch(dataset,
                        all = topALL.mapped,
                        MC = topMC.mapped,
                        HC = topHC.mapped,
                        AC = topAC.mapped,
                        CH= topCH.mapped)
    gene.table <- switch(assembly,
                         GRCh37 = Dixon.ICP.genelist.GRCh37,
                         GRCh38 = Dixon.ICP.genelist.GRCh38)
    
    chrom <- paste0("chr", gene.table$chromosome_name[gene.table$gene == gene])
    start <- gene.table$range_start[gene.table$gene == gene]
    end <- gene.table$range_end[gene.table$gene == gene]
    
    # browser()
    SNP.table$SNP[SNP.table$chr == chrom & SNP.table$bp >= start & SNP.table$bp <= end]
  }
  
  lookup.df <- expand.grid(
    gene = genelist$gene,
    assembly = c("GRCh37"),
    dataset = c("all", "MC", "HC", "AC", "CH")
  )
  lookup.df$search.index <- paste0("s",1:nrow(lookup.df))
  
  # browser()
  overview.list <- apply(lookup.df,1,find_SNP_in_gene_range)
  # browser()
  names(overview.list) <- lookup.df$search.index
  
  overview.list <- overview.list[sapply(overview.list,length) >0]
  
  #swap list keys and values
  
  SNP.gene.matches <- list()
  for (si in names(overview.list)){
    for (j in 1:length(overview.list[[si]])){
      SNP <- as.character(overview.list[[si]][j])
      SNP.gene.matches[[SNP]] <- si
    }
  }
  
  #build lookup lists to annotate SNP.gene.matches with
  # browser()
  top.lookup <- bind_rows(all=topALL.mapped, AC = topAC.mapped, MC = topMC.mapped, HC=topHC.mapped, CH=topCH.mapped, .id = "dataset")
  gene.length.lookup <- bind_rows(GRCh37 = Dixon.ICP.genelist.GRCh37,
                                  # GRCh38 = Dixon.ICP.genelist.GRCh38,
                                  .id = "assembly") %>% mutate(range.length = range_end - range_start) %>% dplyr::select(gene, range.length, assembly)
  
  
  # browser()
  SNP.gene.matches %>% do.call(rbind,.) %>% as.data.frame(.) %>% mutate(SNP = rownames(.)) %>% dplyr::rename("search.index" = V1) %>% 
    left_join(lookup.df, by = "search.index") %>% 
    left_join(top.lookup, by = c("SNP", "dataset")) %>% 
    left_join(gene.length.lookup, by = c("gene", "assembly")) -> SNP.gene.matches
  
  
  BM_anno <- getBM(ensSNP,
                   filters = c("snp_filter"),
                   values = list(SNP.gene.matches$RefSNP),
                   attributes = c("refsnp_id","chrom_start", "refsnp_source","minor_allele_freq"))
  
  SNP.consequences <- getBM(ensSNP,
                            filters = c("snp_filter"),
                            values = list(SNP.gene.matches$RefSNP),
                            attributes = c("refsnp_id","chrom_start", "refsnp_source","minor_allele_freq","consequence_type_tv", "sift_prediction", "sift_score", "polyphen_prediction", "polyphen_score"))
  
  SNP.gene.matches.anno <-merge(SNP.gene.matches, BM_anno, by.x="RefSNP", by.y="refsnp_id", all.x=TRUE)
  SNP.consequences.anno <- merge(SNP.consequences, SNP.gene.matches%>%dplyr::select(RefSNP,gene,chr,bp,OR,`p-value`,rank), by.y="RefSNP", by.x="refsnp_id", all.x=TRUE)
  
  
  
  
  list(SNP.gene.matches=SNP.gene.matches.anno, lookup.df=lookup.df, SNP.consequences=SNP.consequences.anno)
}
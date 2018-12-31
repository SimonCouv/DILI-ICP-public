library(rtracklayer)
library(readr)
library(dplyr)
library(RCurl)
library(XML)
library(stringr)

### get observed mark data for each mark and for both E066 and E118 -------------------------------------------------------------------------------

#load Dixon ranges
Dixon.ICP.genelist.GRCh37 <- read_tsv("Dixon.ICP.genelist.GRCh37.tsv")
Dixon.ICP.genelist.GRCh37 <- Dixon.ICP.genelist.GRCh37 %>% mutate(strand2 = ifelse(Dixon.ICP.genelist.GRCh37$strand > 0,"+","-"),
                                                                  chr=paste0("chr",.$chromosome_name))
#convert to Granges to work with rtracklayer
Dixon.ICP.gr <-Dixon.ICP.genelist.GRCh37 %>% 
  dplyr::select(chr,range_start,range_end,strand2,gene_incl_2kb_prom=gene) %>% 
  makeGRangesFromDataFrame(., 
                           seqnames.field = "chr",
                           start.field = "range_start",
                           end.field = "range_end",
                           keep.extra.columns = TRUE,
                           strand.field = "strand2",
                           seqinfo = Seqinfo(seqnames=NULL, genome = "hg19"))

#convert to BigWigSelection object for use with rtracklayer::import.bw
Dixon.ICP.bw_sel <- BigWigSelection(Dixon.ICP.gr)




### get observed mark data for each mark and for both E066 and E118 -------------------------------------------------------------------------------

get_observed_data <- function(){
  
  #initialize Roadmap Epigenomics data list
  RE.data <- list(
    E066 = list(
      obs = list(),
      imp = list()),
    E118 = list(
      obs= list(),
      imp = list())
  )
  
  print("getting observed data")
  # MACS2 signals most from wustls
  obs.MACS2.url <- "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/"
  marks <- RCurl::getURL(obs.MACS2.url,verbose=FALSE,ftp.use.epsv=TRUE, dirlistonly = TRUE) %>% getHTMLLinks(.) %>% .[grep(.,pattern="^\\w")]
  
  E066.obs.marks <- marks %>% .[grepl(.,pattern = "^E066")] %>% .[grepl(.,pattern = "\\.bigwig$")] %>% .[!grepl(.,pattern = "hg38")]
  E066.obs.marks.links <- as.list(file.path(obs.MACS2.url, E066.obs.marks))
  names(E066.obs.marks.links) <- E066.obs.marks %>% str_match(.,pattern = "-(.*).pval") %>%.[,2]
  
  E118.obs.marks <- marks %>% .[grepl(.,pattern = "^E118")] %>% .[grepl(.,pattern = "\\.bigwig$")] %>% .[!grepl(.,pattern = "hg38")]
  E118.obs.marks.links <- as.list(file.path(obs.MACS2.url, E118.obs.marks))
  names(E118.obs.marks.links) <- E118.obs.marks %>% str_match(.,pattern = "-(.*).pval") %>%.[,2]
  
  
  for (mark in names(E066.obs.marks.links)){
    mark.url <- E066.obs.marks.links[[mark]]
    
    cat("\n")
    print(paste("E066:",mark))
    print(Sys.time())
    
    signal.gr <- import.bw(con = mark.url, selection = Dixon.ICP.bw_sel)
    print(length(signal.gr))
    RE.data[["E066"]][["obs"]][[mark]] <- signal.gr
  }
  
  for (mark in names(E118.obs.marks.links)){
    mark.url <- E118.obs.marks.links[[mark]]
    
    cat("\n")
    print(paste("E118:",mark))
    print(Sys.time())
    
    signal.gr <- import.bw(con = mark.url, selection = Dixon.ICP.bw_sel)
    print(length(signal.gr))
    RE.data[["E118"]][["obs"]][[mark]] <- signal.gr
  }
  
  
  # RNAseq and BSseq signals from wustls
  # https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/README
  # https://egg2.wustl.edu/roadmap/data/byDataType/rna/README
  obs.WGBS.E066 <- "https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/E066_WGBS_FractionalMethylation.bigwig"
  obs.RNAseq.E066p <- "https://egg2.wustl.edu/roadmap/data/byDataType/rna/signal/normalized_bigwig/stranded/E066.Adult_Liver.norm.pos.bw"
  obs.RNAseq.E066n <- "https://egg2.wustl.edu/roadmap/data/byDataType/rna/signal/normalized_bigwig/stranded/E066.Adult_Liver.norm.neg.bw"
  obs.RNAseq.E118p <- "https://egg2.wustl.edu/roadmap/data/byDataType/rna/signal/normalized_bigwig/stranded/E118.HEPG2.norm.pos.bw"
  obs.RNAseq.E118n <- "https://egg2.wustl.edu/roadmap/data/byDataType/rna/signal/normalized_bigwig/stranded/E118.HEPG2.norm.neg.bw"
  
  RE.data[["E066"]][["obs"]][["WGBS"]] <- import.bw(con = obs.WGBS.E066, selection = Dixon.ICP.bw_sel)
  RE.data[["E066"]][["obs"]][["RNAseq_p"]] <- import.bw(con = obs.RNAseq.E066p, selection = Dixon.ICP.bw_sel)
  RE.data[["E066"]][["obs"]][["RNAseq_n"]] <- import.bw(con = obs.RNAseq.E066n, selection = Dixon.ICP.bw_sel)
  RE.data[["E118"]][["obs"]][["RNAseq_p"]] <- import.bw(con = obs.RNAseq.E118p, selection = Dixon.ICP.bw_sel)
  RE.data[["E118"]][["obs"]][["RNAseq_n"]] <- import.bw(con = obs.RNAseq.E118n, selection = Dixon.ICP.bw_sel)
  
  return(RE.data)
}



### get imputed mark data for each mark and for both E066 and E118 -------------------------------------------------------------------------------

get_imputed_data <- function(){
  
  
  #initialize Roadmap Epigenomics data list
  RE.data <- list(
    E066 = list(
      obs = list(),
      imp = list()),
    E118 = list(
      obs= list(),
      imp = list())
  )
  
  RDS_fp <- "../data/rosalind_output/RoadmapEpigenomics_LiverSamples_ICPranges_AllMarks_imp.RDS"
  

  
  print("getting remaining imputed data")
  #remark: RNAseq signal is logRPKM, non-stranded (?)
  
  # get urls to folders listing samples for imputed marks
  imp.url <- "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/"
  imp.marks <- RCurl::getURL(imp.url,verbose=FALSE,ftp.use.epsv=TRUE, dirlistonly = TRUE) %>% getHTMLLinks(.) %>% .[grep(.,pattern="^\\w")] %>% sapply(.,function(x)substr(x,1,nchar(x)-1)[[1]]) %>% as.vector(.)
  imp.mark.links <- as.list(paste0(file.path(imp.url, imp.marks),"/"))
  names(imp.mark.links) <- imp.marks
  
  # loop over marks and their folders
  for (mark in imp.marks){
    

    
    
    if (file.exists(RDS_fp)){
      RE.data <- readRDS(RDS_fp)
      if (mark %in% names(RE.data$E066$imp)|mark %in% names(RE.data$E188$imp)){
        print(sprintf("skipping previously downloaded mark %s", mark))
        next
      } 
    }
    
    cat("\n")
    print(mark)
    print(Sys.time())
    
    #find E066 and E118 files in mark folder
    mark.link <- imp.mark.links[[mark]]
    mark.samples <- RCurl::getURL(mark.link,verbose=FALSE,ftp.use.epsv=TRUE, dirlistonly = TRUE) %>% getHTMLLinks(.) %>% .[grep(.,pattern="^\\w")]
    E066.url <- mark.samples %>% .[grepl(.,pattern = "^E066")] %>% .[grepl(.,pattern = "\\.bigwig$")] %>% .[!grepl(.,pattern = "hg38")] %>% file.path(mark.link,.)
    E118.url <- mark.samples %>% .[grepl(.,pattern = "^E118")] %>% .[grepl(.,pattern = "\\.bigwig$")] %>% .[!grepl(.,pattern = "hg38")] %>% file.path(mark.link,.)
    
    #checks
    stopifnot(length(E066.url) == 1)
    stopifnot(length(E118.url) == 1)
    
    #import data
    signalE066<- import.bw(con = E066.url, selection = Dixon.ICP.bw_sel)
    signalE188<- import.bw(con = E118.url, selection = Dixon.ICP.bw_sel)
    print(c(length(signalE066), length(signalE188)))
    RE.data[["E066"]][["imp"]][[mark]] <- signalE066
    RE.data[["E118"]][["imp"]][[mark]] <- signalE188
    
    saveRDS(object=RE.data, file=RDS_fp)
  }
  
  return(RE.data)
}

RE.data_obs <- get_observed_data()
RE.data_imp <- get_imputed_data()

saveRDS(object=RE.data_obs, file="../data/rosalind_output/RoadmapEpigenomics_LiverSamples_ICPranges_AllMarks_obs.RDS")
saveRDS(object=RE.data_imp, file="../data/rosalind_output/RoadmapEpigenomics_LiverSamples_ICPranges_AllMarks_imp.RDS")






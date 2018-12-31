## ----setup---------------------------------------------------------------
source("https://bioconductor.org/biocLite.R")
source("save_versioned_data.R")
source("read_most_recent.R")



library(dplyr)
library(ggplot2)
library(reshape2)
library(R.utils)
library(tidyr)
library(readr)
library(stringr)
library(data.table)
library(tictoc)


datadir<- file.path(dirname(getwd()),"data")
dbSNPdatadir <- file.path(datadir, "dbSNP_data")
figdir <- file.path(dirname(getwd()),"figs")
RDSdir <- file.path(getwd(),"RDS")

## ----load-data-----------------------------------------------------------
DILI.sum.dir <- file.path(datadir,"Nicoletti_DILI_GWAS_summary_data") 


topAC <- read.table(file.path(DILI.sum.dir,"AC.top.10000.or.txt"), header = FALSE, stringsAsFactors = FALSE) #amoxicillin-clavulanic acid-induced cases
topALL <- read.table(file.path(DILI.sum.dir,"ALL.top.10000.or.txt"), header = FALSE, stringsAsFactors = FALSE) #all patients
topHC <- read.table(file.path(DILI.sum.dir,"Hep.top.10000.or.txt"), header = FALSE, stringsAsFactors = FALSE) #hepatocellular cases
topMC <- read.table(file.path(DILI.sum.dir,"MC.top.10000.or.txt"), header = FALSE, stringsAsFactors = FALSE) #mixed and cholestatic cases
topCH <- read.table(file.path(DILI.sum.dir,"CHOL.10000.top.txt"), header = FALSE, stringsAsFactors = FALSE) #mixed and cholestatic cases

names(topAC) <- names(topALL) <- names(topHC) <- names(topMC) <- names(topCH) <- c("chr", "SNP", "bp", "OR", "95.lower", "95.upper", "p-value")
topAC$rank <- topALL$rank <- topMC$rank <- topHC$rank <- topCH$rank <- 1:10000


## ----data-preparation----------------------------------------------------

Nicoletti.loci <- unique(c(topALL$SNP, topAC$SNP, topHC$SNP, topCH$SNP))
Nicoletti.loci.df <- unique(bind_rows(topALL[,1:3], topAC[,1:3], topHC[,1:3], topMC[,1:3], topCH[,1:3])) %>% mutate(chr_long = paste0("chr", chr))
Nicoletti.loci.DT <- as.data.table(Nicoletti.loci.df)
Nicoletti.loci.DT$chr <- as.character(Nicoletti.loci.DT$chr)
setkey(Nicoletti.loci.DT, chr, bp)
# Nicoletti.loci.df <- union(topALL[,1:3], topAC[,1:3], topHC[,1:3], topMC[,1:3])   #does not work because dplyr::union takes only two arguments


## ----download-and-load-BED-----------------------------------------------
#GRCh37 b151

### write list of ftp urls to download with chono plugin for chrome
write_dbSNP_BED_FTP_URLs <- function(ftp.url,prefix = ""){
  ftp.filenames<- na.omit(str_extract(strsplit(RCurl::getURL(ftp.url), "\r\n")[[1]], pattern = "bed_chr_[\\dXYMT]{1,2}\\.bed.gz$"))
  ftp.filepaths <- file.path(ftp.url,ftp.filenames)
  write.table(ftp.filepaths, paste0(prefix,"ftp.dls.txt"), row.names = FALSE, quote = FALSE)
}


write_dbSNP_BED_FTP_URLs("ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/BED/", prefix = "dbSNPb151_GRCh37p13")
write_dbSNP_BED_FTP_URLs("ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/BED/", prefix = "dbSNPb150_GRCh37p13")
write_dbSNP_BED_FTP_URLs("ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/BED/", prefix = "dbSNPb151_GRCh38p7")
write_dbSNP_BED_FTP_URLs("ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/BED/", prefix = "dbSNPb150_GRCh38p7")

### build list of paths to loop over

humanchromosomes <- c(1:22,"MT","X","Y") %>% .[order(.)]
BEDdir37b151 <- file.path(dbSNPdatadir,"human_9606_b151_GRCh37p13/BED/")
BEDdir37b150 <- file.path(dbSNPdatadir,"human_9606_b150_GRCh37p13/BED/")
BEDdir38b151 <- file.path(dbSNPdatadir,"human_9606_b151_GRCh38p7/BED/")
BEDdir38b150 <- file.path(dbSNPdatadir,"human_9606_b150_GRCh38p7/BED/")
BEDfilespaths37b151 <- file.path(BEDdir37b151,dir(BEDdir37b151))
BEDfilespaths37b150 <- file.path(BEDdir37b150,dir(BEDdir37b150))
BEDfilespaths38b151 <- file.path(BEDdir38b151,dir(BEDdir38b151))
BEDfilespaths38b150 <- file.path(BEDdir38b150,dir(BEDdir38b150))
names(BEDfilespaths37b151) <- names(BEDfilespaths37b150) <- names(BEDfilespaths38b151) <- names(BEDfilespaths38b150) <- humanchromosomes

BEDfilepaths.per.ref <- list(
  GRCh37p13b150 = BEDfilespaths37b150,
  GRCh37p13b151 = BEDfilespaths37b151,
  GRCh38p7b150 = BEDfilespaths38b150,
  GRCh38p7b151 = BEDfilespaths38b151
)


## ----DILI-exploration----------------------------------------------------
#Nicoletti
# table(topALL$chr, useNA = "always")
# table(topAC$chr, useNA = "always")
# table(topHC$chr, useNA = "always")
# table(topMC$chr, useNA = "always")
# table(topCH$chr, useNA = "always")
# 
# #   => Enrichment over chromosomes: clearly more chr6 (HLA) in AC and MC
# 
# venn(list( AC = topAC$SNP, all.pheno = topALL$SNP, CH = topCH$SNP, HC = topHC$SNP))


## ----BED-mapping-FOR-refs-and-codingoffsets------------------------------

rs_map_BED <- function(bedfilepaths.per.ref, codingoffsets){
  gc()
  nicoletti.mapped.all.refs <- list()
  ptm <- proc.time()
  for (ref in names(bedfilepaths.per.ref)){
    # ref <- "GRCh37p13b151" #for debug
    cat(sprintf("started processing %s, elapsed: %s minutes\n", ref,floor((proc.time() - ptm)[3]/60)))
    bedfilepaths <- bedfilepaths.per.ref[[ref]]
    nicoletti.mapped.all.chroms <- list()
    for (mychr in names(bedfilepaths)){
      # mychr <- 22  #for debug
      cat(sprintf("ref %s. started reading chrom %s, elapsed: %s minutes\n", ref, mychr,floor((proc.time() - ptm)[3]/60)))
      chr.BED <- fread(bedfilepaths[mychr], col.names = c("chr", "start", "end", "RefSNP"), 
                       select=1:4, header = FALSE, skip = 1, key = "chr,start,end")
      setkey(chr.BED,chr,start,end)
      nicoletti.mapped.all.offsets <- list()
      for (codingoffset in codingoffsets){
        # codingoffset <- 1 #for debug
        # browser()
        chr.BED[, `:=`(start_w_off = start + codingoffset, end_w_off = end + codingoffset)]
        cat(sprintf("ref %s. starting chrom %s, offset %s. elapsed: %s minutes\n", ref, mychr, codingoffset, floor((proc.time() - ptm)[3]/60)))
        tic()
        nicoletti.mapped.all.offsets[[as.character(codingoffset)]] <- as.data.frame(
          chr.BED[Nicoletti.loci.DT[.(as.character(mychr))], on = .(start_w_off <= bp, end_w_off > bp), .(chr, x.RefSNP, x.start, x.end, bp)]
          # see https://stackoverflow.com/questions/48339895/r-data-table-join-how-to-keep-both-columns-that-are-joined-on
        )
        cat(sprintf("data.table non-equi join:\n"))
        toc()
      }
      nicoletti.mapped.all.chroms[[mychr]] <- bind_rows(nicoletti.mapped.all.offsets, .id = "codingoffset")
      
      rm(chr.BED)
      gc()
    }
    nicoletti.mapped.all.refs[[ref]] <- bind_rows(nicoletti.mapped.all.chroms, .id="chr")
    
  }
  bind_rows(nicoletti.mapped.all.refs, .id = "ref")
}


nicoletti_BED_mapped_full <- rs_map_BED(bedfilepaths = BEDfilepaths.per.ref, codingoffset = (-3):3)
saveRDS(nicoletti_BED_mapped_full, "nicoletti_BED_mapped_full.RDS")

# nicoletti_BED_mapped_direct <- rs_map_BED(bedfilepaths = BEDfilepaths.per.ref[2], codingoffset = 1)
# saveRDS(nicoletti_BED_mapped_direct, "nicoletti_BED_mapped_direct.RDS")







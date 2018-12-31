
#setup -----------------------------------------------------------------------------------------------------------
source("save_versioned_data.R")
source("read_most_recent.R")

library(stringr)
library(biomaRt)
library(dplyr)
#library(profvis)
#library(microbenchmark)
library(ggplot2)
library(reshape2)

library(R.utils)
library(tidyr)
library(readr)
library(rentrez)
library(DBI)
library(dbplyr)
library(data.table)
library(fuzzyjoin)
library(tictoc)

datadir<- file.path(dirname(getwd()),"data")
RDSdir <- file.path(getwd(),"RDS")
rosalind_out_dir <- file.path(datadir,"rosalind_output")
dbSNPdatadir <- file.path(datadir, "dbSNP_data")


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

#load data -----------------------------------------------------------------------------------------------------------

ICP_GWAS_fp <- file.path(datadir,"ICP_GWAS/ICP_notEURonly_log_BRIDGE_PC5.assoc.logistic_send.bz2") 
# ICP_GWAS_DT <- fread(ICP_GWAS_fp,nrows=100,verbose=TRUE) %>% separate(SNP, into = c("chr","bp","allele_A","allele_B"))
cat(sprintf("started reading ICP GWAS data at %s\n", Sys.time()))


ICP_GWAS_DT <- NULL
#encapsulate to deal with idiosyncratic out-of-memory errors on Rosalind, even when allocating sufficient memory on job submission
while (is.null(ICP_GWAS_DT)){
  try(silent=TRUE, expr=
        ICP_GWAS_DT <- fread(ICP_GWAS_fp,verbose=TRUE) %>% separate(SNP, into = c("chr","bp","allele_A","allele_B"))
  )
}
ICP_GWAS_DT[, bp := as.numeric(bp)]

cat(sprintf("Finished pre-processing ICP summary statistics at %s\n", Sys.time()))





#mapping-----------------------------------------------------------------------------------------------------------

ICP_mapped <- list()
ptm <- proc.time()
for (build in c("b150", "b151")){
  for (ref in c("GRCh37p13", "GRCh38p7")){
    
    #debug
    # build= "b150"
    # ref="GRCh37p13"
    
    
    dbSNP_table <- sprintf("rs_%s_%s", build, ref)
    
    for (chrom in 1:22){
      
      #debug
      # chrom=22
      
      bedfilepath <- BEDfilepaths.per.ref[[paste0(ref,build)]][as.character(chrom)]
      chr.BED <- NULL
      cat(sprintf("Started read dbSNP  build %s, ref %s, chrom %s. Wall clock = %s\n", build, ref, chrom, Sys.time()))
      #encapsulate to deal with idiosyncratic out-of-memory errors on Rosalind, even when allocating sufficient memory on job submission
      
      while (is.null(chr.BED)){
        try(silent=TRUE, expr=
              chr.BED <- fread(bedfilepath, col.names = c("chr", "start", "end", "RefSNP"), verbose = TRUE,
                               select=1:4, header = FALSE, skip = 1, key = "chr,start,end"))
      }
      chr.BED <- as.data.table(chr.BED)
      chr.BED[, start:=as.numeric(start)]
      chr.BED[, end:=as.numeric(end)]
    
      for (codingoffset in -3:3){
        
        #debug
        # codingoffset = 1
        
        
        
        
        
        run <- sprintf("%s_%s_%s_%s", build, ref, chrom, sub(as.character(codingoffset),pattern = "-",replacement = "m"))
        
        # run_mapping <-  dbGetQuery(conn = pg_con,
        #                         statement = sprintf("SELECT * FROM icp_gwas LEFT JOIN %s ON (icp_gwas.chr = %s.chr AND %s.startpos + %s <= icp_gwas.bp AND %s.endpos + %s > icp_gwas.bp)",
        #                                             dbSNP_table, dbSNP_table, dbSNP_table, codingoffset, dbSNP_table, codingoffset))
        

        

        
        #print(chrom)

        chr.BED[, `:=`(start_w_off = start + codingoffset, end_w_off = end + codingoffset)]
        
        setkey(chr.BED, RefSNP)
        
        cat(sprintf("build %s ref %s. starting chrom %s, offset %s. elapsed: %s minutes\n", build, ref, chrom, codingoffset, floor((proc.time() - ptm)[3]/60)))
        tic()
        mapped <- chr.BED[ICP_GWAS_DT[chr==(as.character(chrom))], on = .(start_w_off <= bp, end_w_off > bp), .(chr, x.RefSNP, x.start, x.end, bp)]
        # see https://stackoverflow.com/questions/48339895/r-data-table-join-how-to-keep-both-columns-that-are-joined-on
        toc()
        ICP_mapped[[run]] <- mapped
      }
    }
  }
}


#save output -----------------------------------------------------------------------------------------------------------

save_versioned_data(ICP_mapped, RDSdir)
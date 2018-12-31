

################################################################################################################################
# setup ----
################################################################################################################################
library(knitr)
library(kableExtra)
library(gridExtra)
library(stringi)
library(biomaRt)
library(dplyr)
#library(profvis)
#library(microbenchmark)
library(ggplot2)
library(reshape2)
library(xlsx)
library(reader)
library(R.utils)
library(tidyr)
library(readr)
library(stringr)
library(gplots)
library(svMisc)
library(RCurl)
library(data.table)
library(tictoc)
library(rentrez)
library(GenomicRanges)
library(rtracklayer)
library(Gviz)
library(ggman)
library(qqman)
library(podkat)

# point to directories ----------------------------------------------------------------------------------------------------------------------------------------
inputdir <- file.path(dirname(getwd()), "input_data")
datadir<- file.path(dirname(getwd()), "data")
figdir <- file.path(dirname(getwd()), "figs")
RDSdir <- file.path(getwd(),"RDS")
rosalind_outdir <- file.path(RDSdir, "rosalind_output")
DILI.sum.dir <- file.path(inputdir,"Nicoletti_DILI_GWAS_summary_data")

# source function definitions ----------------------------------------------------------------------------------------------------------------------------------------
source("save_versioned_data.R")
source("save_versioned_tsv.R")
source("read_most_recent.R")
source("get_versioned_png_dir.R")  #keeping for compatibility
source("get_versioned_fig_dir.R")
source("inrich_log_parser.R")



# setup BiomaRt -------------------------------------------------------------------------------------------------------------------------
# SNP
ensSNP.GRCh38.p12.mart <- useMart("ENSEMBL_MART_SNP", verbose = TRUE)                                              #ensembl v94
ensSNP.GRCh37.p13.mart <- useMart(biomart = "ENSEMBL_MART_SNP", verbose = TRUE, host = "http://grch37.ensembl.org")                       
ensSNP.GRCh38.p12<- useDataset("hsapiens_snp", mart = ensSNP.GRCh38.p12.mart, verbose = TRUE)
ensSNP.GRCh37.p13<- useDataset("hsapiens_snp", mart = ensSNP.GRCh37.p13.mart, verbose = TRUE)

#genes
ensgenes.GRCh38.p12 <- useMart("ENSEMBL_MART_ENSEMBL" , verbose = TRUE, dataset = "hsapiens_gene_ensembl")
ensgenes.GRCh37.p13 <- useMart("ENSEMBL_MART_ENSEMBL" , verbose = TRUE, dataset = "hsapiens_gene_ensembl", host = "feb2014.archive.ensembl.org")



################################################################################################################################
# load input data ----
################################################################################################################################

#* DILI data ----
topAC <- read.table(file.path(DILI.sum.dir,"AC.top.10000.or.txt"), header = FALSE, stringsAsFactors = FALSE) #amoxicillin-clavulanic acid-induced cases
topALL <- read.table(file.path(DILI.sum.dir,"ALL.top.10000.or.txt"), header = FALSE, stringsAsFactors = FALSE) #all patients
topHC <- read.table(file.path(DILI.sum.dir,"Hep.top.10000.or.txt"), header = FALSE, stringsAsFactors = FALSE) #hepatocellular cases
topMC <- read.table(file.path(DILI.sum.dir,"MC.top.10000.or.txt"), header = FALSE, stringsAsFactors = FALSE) #mixed and cholestatic cases
topCH <- read.table(file.path(DILI.sum.dir,"CHOL.10000.top.txt"), header = FALSE, stringsAsFactors = FALSE) #mixed and cholestatic cases

#Nicoletti et al. 2017
Nicoletti.paper <- data.frame(
  RefSNP = c("rs114577328", "rs72631567", "rs28521457", "rs116561224"),
  overall = c(TRUE, TRUE, FALSE, FALSE),
  MC = c(TRUE, TRUE, FALSE, FALSE),
  HC = c(FALSE, TRUE,TRUE,FALSE),
  drugclass = c(rep("overall",3), "statin")
)

#* ICP data ----
Dixon2014 <- read.csv2(file.path(inputdir, "Dixon2014_AMJG_table1.txt"), header = TRUE)
Dixon.newhits <- readxl::read_xlsx(file.path(inputdir, "Dixon_ICP_dbSNP.xlsx"))
Dixon.ICP.genelist <- readxl::read_xlsx(file.path(inputdir, "Dixon_ICP_candidategenes.xlsx"))

#* annotation data ----
HGNC.table <- readxl::read_xlsx(file.path(datadir, "HGNC_non_alt_loci_set_20181008.xlsx"),col_types = "text")

################################################################################################################################
# process input data ----
################################################################################################################################

#* top-ranking DILI SNPs -------------------------------------------------------------------------------------------------------

#** annotation ----
names(topAC) <- names(topALL) <- names(topHC) <- names(topMC) <- names(topCH) <- c("chr", "SNP", "bp", "OR", "95.lower", "95.upper", "p-value")
topAC$rank <- topALL$rank <- topMC$rank <- topHC$rank <- topCH$rank <- 1:10000

#extend with pure cholestatic hits received 29/10
Nicoletti.loci.ext <- unique(c(topALL$SNP, topAC$SNP, topHC$SNP, topMC$SNP, topCH$SNP))
Nicoletti.loci.ext.df <- unique(bind_rows(topALL[,1:3], topAC[,1:3], topHC[,1:3], topMC[,1:3], topCH[,1:3])) %>% mutate(chr_long = paste0("chr", chr))
Nicoletti.loci.ext.DT <- as.data.table(Nicoletti.loci.ext.df)
setkey(Nicoletti.loci.ext.DT, chr, bp)

#** map top-ranking DILI SNPs to RefSNP ----
nicoletti_BED_mapped_full <- rs_map_BED(bedfilepaths = BEDfilepaths.per.ref, codingoffset = (-3):3)



#* ICP data ----

#** ICP top SNPs ----

#combine RefSNPs from Dixon2014 and newhits WGS
Dixon.RefSNP.data <- bind_rows(Dixon2014 = Dixon2014, Dixon.newhits=Dixon.newhits, .id = "source") %>% arrange(gene, RefSNP)

#see which RefSNPs are in 1KG 
#TODO find query to more specifically filter only for 1KG <<<PHASE 3>>>
Dixon.RefSNP.data$in1KG <- in1KG(Dixon.RefSNP.data$RefSNP)

anno <- getBM(ensSNP.GRCh37.p13,
      filters = "snp_filter",
      values=Dixon.RefSNP.data$RefSNP,
      attributes = c("refsnp_id","chr_name", "chrom_start","chrom_end", "minor_allele_freq")) %>% mutate(RefSNP=refsnp_id)

Dixon.RefSNP.data.GRCh37 <- Dixon.RefSNP.data %>% left_join(anno)


#** ICP genes ----

#join genes listed in the RefSNP table into the candidate gene table
Dixon.ICP.genelist <- full_join(Dixon.ICP.genelist, Dixon.RefSNP.data %>% dplyr::select(gene) %>% distinct())
Dixon.ICP.genelist <- HGNC_annotate(genelist = Dixon.ICP.genelist)


#get gene annotations for Dixon.ICP.genelist, using biomart

################################################################################################################################
# export for external processing ----
################################################################################################################################



################################################################################################################################
# save for data analysis ----
################################################################################################################################

save_versioned_data(Dixon.RefSNP.data.GRCh37, RDSdir)

saveRDS(nicoletti_BED_mapped_full, "nicoletti_BED_mapped_full.RDS")






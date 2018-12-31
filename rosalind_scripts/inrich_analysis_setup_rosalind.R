library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(readr)
library(biomaRt)

datadir<- file.path(dirname(getwd()),"data")
RDSdir <- file.path(getwd(),"RDS")
rosalind_out_dir <- file.path(datadir,"rosalind_output")
rosalind_in_dir <- file.path(RDSdir, "rosalind_output")


source("save_versioned_data.R")
source("read_most_recent.R")

Dixon.ICP.genelist.GRCh37 <- read_tsv("Dixon.ICP.genelist.GRCh37.tsv")


parse_verbose_clump_report <- function(clump_fp){

  clump_con <- file(clump_fp, open = "r")
  clumps <- list()
  clump_membership <- list()
  not_in_reference <- c()
  near_end <- FALSE

  while(TRUE){
    outerline <- str_trim(readLines(clump_con, n=1))
    if(length(outerline)==0 & near_end)
      break

    if (grepl(outerline,pattern="^CHR")){
      clump_active <- TRUE
      clump_snps <- c()
      clump_cnt <-1

      indexline <- readLines(clump_con,n=1)
      index_info <- strsplit(indexline, split=" ")[[1]] %>% .[nchar(.)>0]
      index_chr <- index_info[1]
      index_snp <- index_info[3]
      index_bp <- as.numeric(index_info[4])
      clump_nsnps <- as.numeric(index_info[6])
      clump_nonsig <- as.numeric(index_info[7])
      clump_s05 <- as.numeric(index_info[8])
      clump_s01 <- as.numeric(index_info[9])
      clump_s001 <- as.numeric(index_info[10])
      clump_s0001 <- as.numeric(index_info[11])

      while(clump_active){
        innerline <- readLines(clump_con,n=1)

        # if (length(grepl(innerline,pattern = "RANGE:")) == 0) browser()
        if(grepl(innerline, pattern="RANGE:")){
          # browser()
          rangeline <- innerline
          clump_active <- FALSE
          clump_chr <- str_match(rangeline, pattern="chr(\\d+):")[,2]
          clump_range <- str_match(rangeline, pattern = ":(\\d+\\.\\.\\d+)")[,2]
          clump_start <- as.numeric(strsplit(clump_range, split="\\.\\.")[[1]][1])
          clump_end <- as.numeric(strsplit(clump_range, split="\\.\\.")[[1]][2])
        } else if(grepl(innerline, pattern="------------")){
          # browser()
          clump_active <- FALSE
          clump_chr <- index_chr
          clump_start <- index_bp
          clump_end <- index_bp
        } else {
          if (grepl(innerline,pattern = "rs")){   #if line cont
            clump_snp <- str_match(innerline, pattern = "(rs\\d+)")[,2]
            clump_snps <- c(clump_snps,clump_snp)
          }
        }
      }


      clumps[[index_snp]] <- data_frame(chr = clump_chr,
                                        start = clump_start,
                                        end = clump_end,
                                        span = clump_end - clump_start + 1,
                                        index_chr = index_chr,
                                        index_bp = index_bp,
                                        clump_nsnps = clump_nsnps,
                                        clump_nonsig = clump_nonsig,
                                        clump_s05 = clump_s05,
                                        clump_s01 = clump_s01,
                                        clump_s001 = clump_s001,
                                        clump_s0001 = clump_s0001)

      if(length(clump_snps)>0)
        clump_membership[[index_snp]] <- data_frame(clump_snp = clump_snps)
      else
        clump_membership[[index_snp]] <- data_frame(clump_snp = index_snp)
    }

    if (grepl(outerline, pattern = "not found in dataset")){
      not_in_reference <- c(not_in_reference, str_match(outerline, pattern = "(rs\\d+)")[,2])
      near_end <- TRUE
    }

  }

  clumps.df <- bind_rows(clumps, .id="index_snp")
  clump_membership.df <- bind_rows(clump_membership, .id="index_snp")


  return(list(clumps=clumps.df, clump_membership=clump_membership.df, not_in_reference = not_in_reference))
  close(clump_con)
}

# clump_fp <- "/mnt/lustre/users/k1893262/DILI-ICP/data/VCF_1KG/VCF_1KG_whole_genome/EUR/topAC_EUR_clumped.clumped"
# clump_fp <- "C:/Users/K1893262/OneDrive - King's College London/KCL/DILI_GWAS_Molokhia/data/INRICH/topAC_EUR_clumped.clumped"


# clump_overview_DT <- data.table()
# for (population in c("EUR", "CEU", "GBR")){
# 
#   #debug
#   # population="EUR"
# 
#   print(population)
# 
#   clumpdir <- file.path("/mnt/lustre/users/k1893262/DILI-ICP/data/VCF_1KG/VCF_1KG_whole_genome", population)
#   clump_files <- dir(clumpdir) %>% str_subset(., pattern = "top\\w{2,3}_\\w{3}_\\d_clumped\\.clumped")
# 
#   for (clump_file in clump_files){
# 
#     clump_fp <- file.path(clumpdir, clump_file)
# 
#     rsq_10 <- str_match(clump_file, pattern = "top\\w{2,3}_\\w{3}_(\\d)_clumped\\.clumped")[,2]
#     dili_pheno <- str_match(clump_file, pattern = "top(\\w{2,3})_\\w{3}_\\d_clumped\\.clumped")[,2]
#     print(c(dili_pheno, rsq_10))
#     clump_output <- parse_verbose_clump_report(clump_fp)
#     clump_DT <- data.table(dili_pheno = dili_pheno,
#                            population = population,
#                            rsq_10 = rsq_10,
#                            clumps = list(clump_output$clumps),
#                            membership = list(clump_output$clump_membership),
#                            not_in_reference = list(clump_output$not_in_reference))
#     clump_overview_DT <- rbindlist(list(clump_overview_DT, clump_DT))
#   }
# }

# save_versioned_data(clump_overview_DT, RDSdir)

read_most_recent("clump_overview_DT", RDSdir)




#export input files for inrich  ------------------------------------------------------------------------------------------------

#test interval file
for (population in c("EUR","CEU", "GBR")){
  for (rsq_10 in c(2,5,8)){
    for (dili_pheno in c("AC","ALL", "CH", "HC", "MC")){
      

      clump_overview_DT$clumps[[which(clump_overview_DT$population==population & clump_overview_DT$rsq_10==rsq_10 & clump_overview_DT$dili_pheno==dili_pheno)]] %>%
        dplyr::select(chr=chr, start_bp = start, end_bp = end) %>%
        write_tsv(file.path("/mnt/lustre/users/k1893262/DILI-ICP/data/INRICH/", sprintf("intervals_%s_%s_%s.tsv", dili_pheno, population, rsq_10)))
    }
  }
}


#SNP map file 
for (population in c("EUR","CEU", "GBR")){
  input_bim <- sprintf("/mnt/lustre/users/k1893262/DILI-ICP/data/VCF_1KG/VCF_1KG_whole_genome/%s/%s_allchroms_tmp.bim", population, population)
  target_map <- sprintf("/mnt/lustre/users/k1893262/DILI-ICP/data/INRICH/%s_allchroms.map", population)
  system(sprintf("cut -f1,4 %s > %s", input_bim, target_map))
}

#target set file 
Dixon.ICP.genelist.GRCh37 %>% dplyr::select(ensembl_gene_id) %>% mutate(set_id = "SET1", setname = "ICP_genes") %>%
  write_tsv(col_names = FALSE, path = "/mnt/lustre/users/k1893262/DILI-ICP/data/INRICH/target_set.tsv")

#reference gene file = 19430 protein coding genes from Ensembl BioMart GRCh37, on chromosomes  1-22
ens_genes <- fread("/mnt/lustre/users/k1893262/DILI-ICP/data/INRICH/ens_gene_protcoding_chr1_22_GRCh37.txt")
#reorder columns
ens_genes %>% dplyr::select(4,2,3,1,5) %>% write_tsv(col_names = FALSE, path = "/mnt/lustre/users/k1893262/DILI-ICP/data/INRICH/ens_gene_protcoding_chr1_22_GRCh37_clean.tsv")













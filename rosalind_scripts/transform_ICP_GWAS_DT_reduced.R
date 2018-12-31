library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
source("save_versioned_data.R")
source("get_versioned_fig_dir.R")

RDSdir <- file.path(getwd(),"RDS")

ICP_GWAS_mapped_DT <- readRDS("/mnt/lustre/users/k1893262/DILI-ICP/DILI-ICP/RDS/ICP_GWAS_mapped_DT_v20181214_132346.RDS")
# ICP_GWAS_mapped_DT <- rbindlist(ICP_GWAS_mapped, idcol = "run")
# ICP_GWAS_mapped_DT[, c("build","ref","chr","codingoffset") := tstrsplit(run, "_", fixed=TRUE)]
# ICP_GWAS_mapped_DT %>% separate(col=run, into=c("build","ref","chr","offset"),sep="_") %>% as.data.table(.) -> ICP_GWAS_mapped_DT_sep
# save_versioned_data(ICP_GWAS_mapped_DT, RDSdir)

p.ICP_GWAS_mapping_yield_chr22 <- ICP_GWAS_mapped_DT[chr=="22"] %>%
  dplyr::rename_all(~sub('x\\.', '', .x)) %>% 
  group_by(chr, bp, ref, build, codingoffset) %>% 
  summarise(n.RefSNP.per.DILI.SNP = sum(!is.na(RefSNP))) %>%  
  group_by(ref,build, codingoffset) %>% 
  summarise(yield=mean(n.RefSNP.per.DILI.SNP>0)) %>%
  
  ggplot(aes(x=interaction(build,ref),y=yield,fill=codingoffset)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw()+
  xlab("Reference genome")+
  ylab("Yield")+
  scale_fill_discrete(name="offset")
save_versioned_data(p.ICP_GWAS_mapping_yield_chr22, RDSdir)

ICP_GWAS_mapped_GRCh37_b151_DT2 <- ICP_GWAS_mapped_DT[ref=="GRCh37p13" & build=="b151"]
# save_versioned_data(ICP_GWAS_mapped_GRCh37_b151_DT)
save_versioned_data(ICP_GWAS_mapped_GRCh37_b151_DT2, RDSdir)
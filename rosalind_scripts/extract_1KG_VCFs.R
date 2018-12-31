#https://www.biostars.org/p/2909/#75824

library(readr)
Dixon.ICP.genelist.GRCh37 <- read_tsv("Dixon.ICP.genelist.GRCh37.tsv")

#system("module load bioinformatics/samtools/1.5") does not work

for (i in 1:nrow(Dixon.ICP.genelist.GRCh37)){
  gene = Dixon.ICP.genelist.GRCh37$gene[i]
  chr = Dixon.ICP.genelist.GRCh37$chromosome_name[i]
  start = Dixon.ICP.genelist.GRCh37$range_start[i]
  end = Dixon.ICP.genelist.GRCh37$range_end[i]
  print(gene)
  system(sprintf(
"tabix -h ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz %s:%s-%s > /users/k1893262/brc_scratch/DILI-ICP/data/VCF_1KG_ICP_regions/full_VCFs/%s_1KG.vcf",
chr,chr,start,end,gene))
}

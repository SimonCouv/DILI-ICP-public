#!/bin/bash

pop="$1"
samplesfile=/users/k1893262/brc_scratch/DILI-ICP/data/VCF_1KG_ICP_regions/sample_populations/"$pop"_samples.tsv
mkdir -p /users/k1893262/brc_scratch/DILI-ICP/data/VCF_1KG_ICP_regions/"$pop"/

for f in /users/k1893262/brc_scratch/DILI-ICP/data/VCF_1KG_ICP_regions/full_VCFs/*.vcf; do 
  gene=$(echo $f | grep -oP "(\w+)(?=_1KG.vcf)")
  echo "$gene"
  #zcat "$f" | grep "^.[^#]" | head -2 | cut -c1-20 > ~/rds/hpc-work/chr"$chr"_extracted.vcf
  #echo $samples
  bcftools view -O z -S "$samplesfile" "$f" > /users/k1893262/brc_scratch/DILI-ICP/data/VCF_1KG_ICP_regions/"$pop"/"$gene"_1KG_"$pop".vcf
  #echo "success"
done



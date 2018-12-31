#!/bin/bash
#$ -pe smp 4
#$ -cwd
#$ -N plink_clump
#$ -j y
#$ -l h_vmem=19G
#$ -l h_rt=24:00:00
#$ -o /mnt/lustre/users/k1893262/DILI-ICP/DILI-ICP/rosalind_logs/
#$ -R y
#$ -M simon.couvreur@kcl.ac.uk
#$ -m be
module load bioinformatics/plink2/1.90b3.38
module load bioinformatics/R/3.5.0 







DIR1KG=/mnt/lustre/users/k1893262/DILI-ICP/data/VCF_1KG/VCF_1KG_whole_genome/ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/
DILI_DIR=/mnt/lustre/users/k1893262/DILI-ICP/data/Nicoletti_DILI_GWAS_summary_data/
olddir=($pwd)


#download all whole genome VCFs
#wget -r ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/

for pop in EUR CEU GBR; do
	samplesfile=/users/k1893262/brc_scratch/DILI-ICP/data/VCF_1KG/sample_populations/"$pop"_samples.tsv
	DIRVCF=/users/k1893262/brc_scratch/DILI-ICP/data/VCF_1KG/VCF_1KG_whole_genome/"$pop"/
	
	mkdir -p $DIRVCF
	touch "$DIRVCF"/bfile_list.txt
	cd $DIRVCF
	
	for CHR in {1..22}; do

	   VCF_IN="$DIR1KG"/ALL.chr"$CHR".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	   

	   # #subset VCFs using plink
	   # #https://adairama.wordpress.com/2017/05/12/download-1000-genomes-phase3-and-calculate-allele-frequencies/	

	   #plink --vcf $VCF_IN --biallelic-only strict --keep-fam $samplesfile --geno 0.05 --hwe 1e-50 --make-bed --out "$pop"_chr"$CHR"_tmp    
			#https://www.cog-genomics.org/plink/1.9/filter
			#biallelic-only only to avoid downstream merge problems (https://www.cog-genomics.org/plink/1.9/input#vcf and https://www.cog-genomics.org/plink/1.9/data#merge3)
	   plink --bfile "$pop"_chr"$CHR"_tmp --exclude "$pop"_allchroms_tmp-merge.missnp --make-bed --out "$pop"_chr"$CHR"_filter
			#should not be needed with --biallelic-only strict flag
	   
	   echo "$pop"_chr"$CHR"_filter.bed "$pop"_chr"$CHR"_filter.bim "$pop"_chr"$CHR"_filter.fam >> bfile_list.txt
	   #echo "$pop"_chr"$CHR"_tmp.bed "$pop"_chr"$CHR"_tmp.bim "$pop"_chr"$CHR"_tmp.fam >> bfile_list.txt
	done



	#merge plink binary filesets
	#plink --bfile "$pop"_chr1_tmp --merge-list bfile_list.txt --make-bed --out "$pop"_allchroms_tmp
		#generates error on first run + misSNP file which can be used with --exclude
	# https://www.biostars.org/p/116960/

	plink --bfile "$pop"_chr1_filter --merge-list bfile_list.txt --make-bed --out "$pop"_allchroms_tmp

	#clumping step
	#http://zzz.bwh.harvard.edu/plink/clump.shtml
	#https://www.cog-genomics.org/plink/1.9/postproc#clump

	for rsq_10 in {1..9};do
		rsq=$(Rscript -e "cat("$rsq_10"/10)")

		ls "$DILI_DIR" | grep ".mapped.assoc$" | while read -r file ; do
			echo "Clumping $file with 1KG $pop as LD reference and $rsq as r2 threshold."
			PHENOTYPE=$(grep -oP '\w+(?=_mapped.assoc)' <<< $file)
			ASSOCPATH="$DILI_DIR""$file"
			REPORT_OUT=/users/k1893262/brc_scratch/DILI-ICP/data/VCF_1KG/VCF_1KG_whole_genome/"$pop"/"$PHENOTYPE"_"$pop"_"$rsq_10"_clumped
			plink --bfile "$pop"_allchroms_tmp --clump $ASSOCPATH --clump-verbose --clump-p1 0.0003 --clump-r2 $rsq --clump-kb 250 --out $REPORT_OUT
		done
	done
done




#rm *_tmp.*
#rm EUR_$CHR.{log,nosex}

cd $olddir


#plink --file mydata --clump mytest1.assoc



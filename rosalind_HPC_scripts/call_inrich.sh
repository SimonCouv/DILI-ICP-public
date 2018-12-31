#!/bin/sh
#$ -pe smp 4
#$ -cwd
#$ -N inrich
#$ -j y
#$ -o ./rosalind_logs/
#$ -l h_vmem=19G
#$ -l h_rt=24:00:00
#$ -R y
#$ -M simon.couvreur@kcl.ac.uk
#$ -m be

module load bioinformatics/R/3.5.0 
PATH=$PATH:/users/k1893262/bin/inrich
export PATH

cd /mnt/lustre/users/k1893262/DILI-ICP/data/INRICH/

for line in {1..45}; do

	INTERVALFILE=$(awk "NR==$line" intervals_list.txt)
	echo $INTERVALFILE
	
	#following line caused code not to work for dili_pheno ALL
	#POP=$(Rscript -e "cat(stringr::str_match('"$INTERVALFILE"', pattern='intervals_\\\\w{2}_(\\\\w{2,3})_\\\\d\\\\.tsv')[,2])")
	
	#fixed:
	POP=$(Rscript -e "cat(stringr::str_match('"$INTERVALFILE"', pattern='intervals_\\\\w{2,3}_(\\\\w{3})_\\\\d\\\\.tsv')[,2])")

	inrich -a "$INTERVALFILE" -m "$POP"_allchroms.map -t target_set.tsv -g ens_gene_protcoding_chr1_22_GRCh37_clean.tsv -p 1 -o inrich_out_"$line"

done
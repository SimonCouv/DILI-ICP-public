#!/bin/bash

module load bioinformatics/plink2/1.90b3.38

pop="$1"
pattern="([\w\d]*)(?=_1KG[\w\d_]*.vcf$)"
regionpath=/users/k1893262/brc_scratch/DILI-ICP/data/VCF_1KG_ICP_regions

genes=$(ls "$regionpath"/"$pop"/ | grep -oP "$pattern")
# pattern="[.]*\.vcf$"
# greprexeg="([\w\d]*)(?=.vcf$)"
for gene in $genes
do
	echo "$gene"
	fp="$regionpath"/"$pop"/"$gene"_1KG_"$pop".vcf
	#ls "$fp"
	plink --vcf "$fp" --double-id --r2 inter-chr gz dprime with-freqs --ld-window-r2 0 --out "$regionpath"/"$pop"/"$gene"_"$pop"

done

oldwd=$(pwd)

cd "$regionpath"/"$pop"/
echo *.ld.gz | xargs tar -cf "$pop"_LD.tar
rm *.ld.gz

# (cd "$regionpath"/"$pop"/ && tar - cvzf *.ld.gz ) | gzip -9 > "$regionpath"/"$pop"/"$pop"_LD.gz
# rm "$regionpath"/"$pop"/*.ld.gz


# for vcf in /mnt/lustre/users/k1893262/DILI-ICP/data/VCF_1KG_ICP_regions/*.vcf; do
    # vcf=
    # out=
    # plink --vcf "$vcf" --double-id --r2 inter-chr dprime with-freqs --out
# done


cd $oldwd

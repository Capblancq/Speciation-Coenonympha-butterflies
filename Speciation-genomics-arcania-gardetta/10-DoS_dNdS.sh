#!/bin/bash

#SBATCH --job-name=FixedDivergence
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=20G
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

## Source required programs
. /local/env/envvcftools-0.1.16.sh
. /local/env/envbcftools-1.9.sh
. /local/env/envperl-5.26.2.sh
. /local/env/envconda.sh
conda activate ~/TOOLS/vcflib

## Output folder
input=/scratch/tcapblancq/Genotypes
output=/scratch/tcapblancq/dNdS
cd ${output}

## The vcf file had to be bgzipped and indexed as below before running that script 
bcftools view -I ${input}/All_polymorphicsites_annotated.vcf -O z -o ${input}/All_polymorphicsites_annotated.vcf.gz
bcftools index ${input}/All_polymorphicsites_annotated.vcf.gz

for i in {1..30}
do

	# Target chromosome
	CHROM="scaffold_${i}"

	# Subset vcf based on species and scaffold i
	bcftools view -s ind0,ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8,ind9,ind10,ind11,ind12,ind13,ind14,ind15,ind16,ind17,ind18 -r ${CHROM} ${input}/All_polymorphicsites_annotated.vcf.gz > Carcania_${CHROM}.vcf
	bcftools view -s ind19,ind20,ind21,ind22,ind23,ind24,ind25,ind26,ind27,ind28,ind29,ind30,ind31,ind32,ind33,ind34 -r ${CHROM} ${input}/All_polymorphicsites_annotated.vcf.gz > Cgardetta_${CHROM}.vcf

	# Re-estimate allele frequencies in the subsetted vcf files and keep only the first 8 columns
 	vcffixup Carcania_${CHROM}.vcf | cut -f1-8 | grep 'synonymous\|missense' > Carcania_coding_${CHROM}.vcf
 	vcffixup Cgardetta_${CHROM}.vcf | cut -f1-8 | grep 'synonymous\|missense' > Cgardetta_coding_${CHROM}.vcf

 	# Remove temporary files
 	rm Carcania_${CHROM}.vcf
 	rm Cgardetta_${CHROM}.vcf

done

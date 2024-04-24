#!/bin/bash

#SBATCH --job-name=coverageBAM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=10G
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

# Load the required programs and packages
. /local/env/envsamtools-1.15.sh

# Working Directories
input=/groups/divalps/Mapping/WGS/bam/
mkdir /scratch/tcapblancq/Sequencing
output=/scratch/tcapblancq/Sequencing

# To do first
awk 'NR>1{print $1,$2-25000,$2+25000}' windows.txt > windows.bed

# Estimates coverage for each window along the genome
for BAM in ${input}/*.final.bam
do
	sample=${BAM/.final.bam/}
	name=`basename ${sample}`
	
	# Find per base coverage for each window
	samtools bedcov windows.bed ${BAM} -Q 20 -d 4 > ${output}/coverage_${name}.txt
done

# Find the number of SNP per window
while IFS= read -r line
do
	# Get the information about the window
	CHROM=$( echo ${line} | awk '{print $1}' )
	POS=$( echo ${line} | awk '{print $2}' )
	
	# Number of SNPs on that window
	NB_SNP_arc=`zcat /scratch/tcapblancq/Genotypes/arcania_polymorphicsites.mafs.gz | awk -v chrom=${CHROM} -v pos=${POS} 'FNR>1{if ($1==chrom && $2>=(pos-25000) && $2<=(pos+25000)) print $0}' | wc -l`
        NB_SNP_gar=`zcat /scratch/tcapblancq/Genotypes/gardetta_polymorphicsites.mafs.gz | awk -v chrom=${CHROM} -v pos=${POS} 'FNR>1{if ($1==chrom && $2>=(pos-25000) && $2<=(pos+25000)) print $0}' | wc -l`

	touch ${output}/NB_snps_windows_arcania.txt
	echo ${CHROM} ${POS} ${NB_SNP_arc} >> ${output}/NB_snps_windows_arcania.txt 

        touch ${output}/NB_snps_windows_gardetta.txt
        echo ${CHROM} ${POS} ${NB_SNP_gar} >> ${output}/NB_snps_windows_gardetta.txt

done < windows2.txt

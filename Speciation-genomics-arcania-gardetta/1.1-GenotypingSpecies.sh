#!/bin/bash

#SBATCH --job-name=SNPcallingANGSD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH --array=0-1
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

# Defining the file to process using the Array ID
files=(~/Study-Speciation-arcania-gardetta/*a_bam.list)
file=${files[$SLURM_ARRAY_TASK_ID]}

# Retreive the name of the species from the file name
species=${file/_bam.list}
name=`basename ${species}`

# Source required programs
. /local/env/envconda.sh
conda activate ~/TOOLS/angsd-0.921

# Output folder
output=/scratch/tcapblancq/Genotypes

# Reference genome
ref=~/Study-Mapping-WGS/genome/BST1.fa

# Estimating genotype likelihoods for the individuals included in the *_bam.list
angsd -b ${file} \
	-ref ${ref} \
	-anc ${ref} \
	-out ${output}/${name}_allsites \
	-nThreads 10 \
	-uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
	-setMinDepthInd 4 -setMaxDepthInd 80 -minInd 5 \
	-skipTriallelic 1 \
	-GL 1 \
	-doMajorMinor 1 -doMaf 1 \
	-doCounts 1

# Finding loci intersection among regions
#zcat ${output}/arcania_allsites.mafs.gz | awk 'BEGIN { OFS = ":" }{ print $1,$2 }' | sed '1d' | sort > ${output}/arcania_allsites_sites.txt
#zcat ${output}/gardetta_allsites.mafs.gz | awk 'BEGIN { OFS = ":" }{ print $1,$2 }' | sed '1d' | sort > ${output}/gardetta_allsites_sites.txt

#comm -12 ${output}/arcania_allsites_sites.txt ${output}/gardetta_allsites_sites.txt > ${output}/arcania_gardetta_allsites_sites.txt

#sed 's/:/\t/' ${output}/arcania_gardetta_allsites_sites.txt | sort -b -k1,1 > ${output}/intersect_arcgar_allsites.txt
#cut -f1 ${output}/intersect_arcgar_allsites.txt | uniq | sort > ${output}/intersect_arcgar_allsites.chrs
#angsd sites index ${output}/intersect_arcgar_allsites.txt

# Remove intermediary files
#rm ${output}/arcania_allsites_sites.txt
#rm ${output}/gardetta_allsites_sites.txt
#rm ${output}/arcania_gardetta_allsites_sites.txt

# Estimating genotype likelihoods for the individuals included in the *_bam.list but only intersecting sites
angsd -b ${file} \
	-ref ${ref} \
	-anc ${ref} \
	-out ${output}/${name}_allsites \
	-sites ${output}/intersect_arcgar_allsites.txt \
	-nThreads 10 \
	-uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
	-GL 1 \
	-doMaf 1 -doMajorMinor 1 -doSaf 1

# Estimating genotype likelihoods for the individuals included in the *_bam.list but only intersecting polymorphic sites
angsd -b ${file} \
        -ref ${ref} \
        -anc ${ref} \
        -out ${output}/${name}_polymorphicsites \
        -sites ${output}/intersect_arcgar_allsites.txt \
        -nThreads 10 \
        -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
        -GL 1 \
        -snp_pval 1e-6 \
	-doMajorMinor 1 -doMaf 1 -doGlf 2

#!/bin/bash

#SBATCH --job-name=GenotypingANGSD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=30G
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

# Source required programs
. /local/env/envconda.sh
conda activate ~/TOOLS/angsd-0.921

# Output folder
output=/scratch/tcapblancq/Genotypes

# Reference genome
ref=~/Study-Mapping-WGS/genome/BST1.fa

# Estimating genotype likelihoods for the individuals included in the *_bam.list
angsd -b ~/Study-Speciation-arcania-gardetta/all_bam.list \
	-ref ${ref} \
	-anc ${ref} \
	-out ${output}/All_polymorphicsites \
        -sites ${output}/intersect_arcgar_allsites.txt \
        -nThreads 10 \
        -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
        -GL 1 \
        -snp_pval 1e-6 \
	-doMaf 1 -doMajorMinor 1 -doGlf 2 -dopost 1 -doGeno 1 -postCutoff 0.8 -doVcf 1

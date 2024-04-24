#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=72:00:00
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

# Output folder
input=/scratch/tcapblancq/Genotypes
output=/scratch/tcapblancq/LD

# Creating a pos file
zcat ${input}/${name}_polymorphicsites.mafs.gz | awk '{print $1"\t"$2}' | sed '1d' > ${output}/pos_${name}.txt

# Estimating LD
count=`cat ${file} | wc -l`
N_SITES=$((`zcat ${input}/${name}_polymorphicsites.mafs.gz | wc -l`-1))
~/TOOLS/ngsLD/ngsLD --geno ${input}/${name}_polymorphicsites.beagle.gz --probs --n_ind ${count} --n_sites ${N_SITES} --pos ${output}/pos_${name}.txt --out ${output}/${name}.LD --min_maf 0.05 --rnd_sample 0.5 --max_kb_dist 10

# Removing temporary files
rm ${output}/pos_${name}.txt

#!/bin/bash

#SBATCH --job-name=DILS_windows
#SBATCH --cpus-per-task=10
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=10G
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

# Output folder
output=/scratch/tcapblancq/DILS

## Concatenate all the file stogether (has to be done once all the window files are done!)
find ${output}/ -size 0 -delete
rm ${output}/input.dils.fasta
cat ${output}/*.dils > ${output}/input.dils.fasta
rm ${output}/*.dils

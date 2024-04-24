#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=30G
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

# Source required programs
. /local/env/envconda.sh
conda activate ~/TOOLS/angsd-0.921

# Output folder
input=/scratch/tcapblancq/Genotypes
output=/scratch/tcapblancq/Thetas

# Estimating thetas from the SFS
realSFS ${input}/arcania_allsites.saf.idx -maxIter 50000 -tole 1e-6 -P 10 > ${output}/arcania.sfs
realSFS saf2theta ${input}/arcania_allsites.saf.idx -outname ${output}/arcania -sfs ${output}/arcania.sfs
thetaStat do_stat ${output}/arcania.thetas.idx -win 50000 -step 50000 -type 0 -outnames ${output}/arcania

realSFS ${input}/gardetta_allsites.saf.idx -maxIter 50000 -tole 1e-6 -P 10 > ${output}/gardetta.sfs
realSFS saf2theta ${input}/gardetta_allsites.saf.idx -outname ${output}/gardetta -sfs ${output}/gardetta.sfs
thetaStat do_stat ${output}/gardetta.thetas.idx -win 50000 -step 50000 -type 0 -outnames ${output}/gardetta


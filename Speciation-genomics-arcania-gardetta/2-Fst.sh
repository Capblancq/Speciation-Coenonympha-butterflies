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
output=/scratch/tcapblancq/Fst

# Optimization of the 2D SFS for the four species
realSFS ${input}/arcania_allsites.saf.idx ${input}/gardetta_allsites.saf.idx -maxIter 50000 -tole 1e-6 -P 10 > ${output}/arcania.gardetta.sfs

# Prepare the Fst
realSFS fst index ${input}/arcania_allsites.saf.idx ${input}/gardetta_allsites.saf.idx \
	-sfs ${output}/arcania.gardetta.sfs \
	-fstout ${output}/arcania-gardetta_Fst

# Get the estimate along sliding windows
realSFS fst stats2 ${output}/arcania-gardetta_Fst.fst.idx -win 50000 -step 50000 -type 0 > ${output}/Fst_sliding_arcania-gardetta.txt

# Get Fst estimate for each site
realSFS fst print ${output}/arcania-gardetta_Fst.fst.idx > ${output}/Fst_sites_arcania-gardetta.txt

#!/bin/bash

#SBATCH --job-name=Fastsimcoal2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --array=1-30%5
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

# Run ID
i=$SLURM_ARRAY_TASK_ID

# Scenario 
PREFIX=AM

# Move to working directory
cd /scratch/tcapblancq/Fastsimcoal

# Create directory for the scenario
mkdir ${PREFIX}
mkdir ${PREFIX}/run${i}

# Copy and rename all required files
cp ${PREFIX}.tpl ${PREFIX}.est ${PREFIX}/run${i}/
cp Coenonympha_jointMAFpop1_0.obs ${PREFIX}/run${i}/${PREFIX}_jointMAFpop1_0.obs

# Enter the run directory
cd ${PREFIX}/run${i}

# Run fastsimcaol
~/TOOLS/fsc2705 -t ${PREFIX}.tpl -n1000000 -m -e ${PREFIX}.est -M -L40 -q -c10

#!/bin/bash

#SBATCH --job-name=Fastsimcoal2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --array=30-100%10
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

# Defining the file to process using the Array ID
i=${SLURM_ARRAY_TASK_ID}
echo ${i}

# Input/output directory
DIR=/scratch/tcapblancq/Fastsimcoal
cd ${DIR}/bootstrap
SCENARIO=SC

# Running the optimization with each bootstraped SFS
cp ${SCENARIO}.tpl ./${SCENARIO}/${SCENARIO}_$i/${SCENARIO}.tpl
cp ${SCENARIO}.est ./${SCENARIO}/${SCENARIO}_$i/${SCENARIO}.est
cp ${SCENARIO}.pv ./${SCENARIO}/${SCENARIO}_$i/${SCENARIO}.pv
cd ./${SCENARIO}/${SCENARIO}_$i/
~/TOOLS/fsc2705 -t ${SCENARIO}.tpl -e ${SCENARIO}.est -n1000000 -d -M -L50 --initValues ${SCENARIO}.pv -c10 -q

#!/bin/bash

#SBATCH -p genscale
#SBATCH --job-name=DILS_inference
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=20G
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

# Source required programs
#. /local/env/envpython-3.7.1.sh
. ~/TOOLS/TUTU/bin/activate
. /local/env/envsnakemake-6.0.5.sh
. /local/env/envperl-5.26.2.sh
. /local/env/envconda.sh
conda activate ~/TOOLS/pypy
. /local/env/envr-4.1.3.sh

# Working directory
DIR=/scratch/tcapblancq/DILS

# DILS binaries
DILS=/home/genouest/cnrs_umr5553/tcapblancq/TOOLS/DILS/bin

# Runs DILS
snakemake --snakefile ${DILS}/Snakefile_2pop -p -j 140 --configfile ${DIR}/myConfigFileTwoPops_noZ.yaml --cluster-config ${DILS}/cluster_2pop.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu=30G" --latency-wait 60 --rerun-incomplete #--unlock


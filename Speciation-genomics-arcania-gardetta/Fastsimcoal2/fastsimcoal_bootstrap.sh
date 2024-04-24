#!/bin/bash

#SBATCH -p genscale
#SBATCH --job-name=Fastsimcoal2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=10G
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

# Input/output directory
DIR=/scratch/tcapblancq/Fastsimcoal

# Creates a sub-directory
#mkdir ${DIR}/bootstrap
cd ${DIR}/bootstrap

# Have to find the .par file of the best run using R treatment of the results
SCENARIO=SC
#RUN=16
#cp ${DIR}/${SCENARIO}/run${RUN}/${SCENARIO}.par ./
#cp ${DIR}/${SCENARIO}/run${RUN}/${SCENARIO}/${SCENARIO}.pv ./
#cp ${DIR}/${SCENARIO}.tpl ./
#cp ${DIR}/${SCENARIO}.est ./

# Have to change a bit the .par file to 
#sed -i 's/^1 0$/200000 0/g' ${SCENARIO}.par 
#sed -i 's/^FREQ 1/DNA 100/g' ${SCENARIO}.par

# Then generate 100 SFS
#~/TOOLS/fsc2705 -i ${SCENARIO}.par -n100 -j -d -s0 -x -I -q -c10

# Running the optimization with each bootstraped SFS
#for i in {17..100}
#do
#	cp ${SCENARIO}.tpl ./${SCENARIO}/${SCENARIO}_$i/${SCENARIO}.tpl
#	cp ${SCENARIO}.est ./${SCENARIO}/${SCENARIO}_$i/${SCENARIO}.est
#	cp ${SCENARIO}.pv ./${SCENARIO}/${SCENARIO}_$i/${SCENARIO}.pv
#	cd ./${SCENARIO}/${SCENARIO}_$i/
#	~/TOOLS/fsc2705 -t ${SCENARIO}.tpl -e ${SCENARIO}.est -n1000000 -d -M -L50 --initValues ${SCENARIO}.pv -c10 -q
#	cd ${DIR}/bootstrap/
#done

# Estimating confidence intervals
cat ./${SCENARIO}/${SCENARIO}_1/${SCENARIO}/${SCENARIO}.bestlhoods | awk 'NR==1' > Parameters_bootstrap.txt
for i in {1..16}
do
	cat ./${SCENARIO}/${SCENARIO}_$i/${SCENARIO}/${SCENARIO}.bestlhoods | awk 'NR==2' >> Parameters_bootstrap.txt
done

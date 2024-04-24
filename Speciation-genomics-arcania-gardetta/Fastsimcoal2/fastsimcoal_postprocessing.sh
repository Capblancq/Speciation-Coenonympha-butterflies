#!/bin/bash

#SBATCH --job-name=Fastsimcoal2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=20G
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null


## Retreiving results from all scenario and iterations
for PREFIX in IM SC AM SI #TM
do 
	echo "run" `cat /scratch/tcapblancq/Fastsimcoal/${PREFIX}/run1/${PREFIX}/${PREFIX}.bestlhoods | awk 'NR==1'` > results_${PREFIX}.txt
	for i in {1..30}
	do
		if test -f /scratch/tcapblancq/Fastsimcoal/${PREFIX}/run${i}/${PREFIX}/${PREFIX}.bestlhoods
		then
			echo "$i" `cat /scratch/tcapblancq/Fastsimcoal/${PREFIX}/run${i}/${PREFIX}/${PREFIX}.bestlhoods | awk 'NR==2'` >> results_${PREFIX}.txt
		fi
	done
done

#!/bin/bash

#SBATCH --job-name=LD_windows
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --array=3001-9843%50
#SBATCH        --error=/dev/null
#SBATCH        --output=/dev/null

##SLURM_ARRAY_TASK_ID=$((${SLURM_ARRAY_TASK_ID}+10000))

# Extract chrom and pos of the target window
CHROM=$( sed -n $((${SLURM_ARRAY_TASK_ID}+1))p windows.txt | awk '{split($1,a,"_")}{print a[2]}' )
POS=$( sed -n $((${SLURM_ARRAY_TASK_ID}+1))p windows.txt | awk '{print $2}' )

# Input/Output folders
INPUT=/scratch/tcapblancq/Genotypes
OUTPUT=/scratch/tcapblancq/LD

# Source R
. /local/env/envconda.sh
conda activate ~/TOOLS/R

for species in arcania gardetta
do

  # Subsetting the beagle files based on the windows
  zcat ${INPUT}/${species}_polymorphicsites.beagle.gz | awk -v chrom=${CHROM} -v pos=${POS} 'FNR==1{print}{split($1,a,"_")}{if (a[2]==chrom && a[3]>=(pos-25000) && a[3]<=(pos+25000)) print $0}' | gzip > ${OUTPUT}/${species}_window${SLURM_ARRAY_TASK_ID}.beagle.gz

  # Creating a pos file
  zcat ${OUTPUT}/${species}_window${SLURM_ARRAY_TASK_ID}.beagle.gz | sed '1d' | awk '{split($1,a,"_")}{print a[1]"_"a[2]"\t"a[3]}' > ${OUTPUT}/pos_${species}_window${SLURM_ARRAY_TASK_ID}.txt

  # Estimating LD
  COUNT=`cat ${species}_bam.list | wc -l`
  N_SITES=$((`cat ${OUTPUT}/pos_${species}_window${SLURM_ARRAY_TASK_ID}.txt | wc -l`))

  # Retreive mean r2 on the window 
  if (("${N_SITES}" > "10"))
  then
	#Estime LD with ngsLD
	~/TOOLS/ngsLD/ngsLD --geno ${OUTPUT}/${species}_window${SLURM_ARRAY_TASK_ID}.beagle.gz --probs --n_ind ${COUNT} --n_sites ${N_SITES} --pos ${OUTPUT}/pos_${species}_window${SLURM_ARRAY_TASK_ID}.txt --out ${OUTPUT}/${species}_window${SLURM_ARRAY_TASK_ID}.LD

	# Estimate r2 per base
	R2=`awk '{ print $4/$3 }' ${OUTPUT}/${species}_window${SLURM_ARRAY_TASK_ID}.LD | awk '{sum+=$1} END { print sum/NR}'`	

	# Feeding the output file
 	echo ${CHROM} ${POS} ${R2} ${N_SITES} >> ${OUTPUT}/${species}_LDr2bp.txt
  else
	echo ${CHROM} ${POS} "NA" ${N_SITES} >> ${OUTPUT}/${species}_LDr2bp.txt
  fi

  # Removing temporary files
  rm ${OUTPUT}/${species}_window${SLURM_ARRAY_TASK_ID}.beagle.gz
  rm ${OUTPUT}/pos_${species}_window${SLURM_ARRAY_TASK_ID}.txt
  rm ${OUTPUT}/${species}_window${SLURM_ARRAY_TASK_ID}.LD

done



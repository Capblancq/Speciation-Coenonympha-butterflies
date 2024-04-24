#!/bin/bash

#SBATCH --job-name=DILS_windows
#SBATCH --cpus-per-task=10
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --array=9001-9844%20
#SBATCH        --error=/dev/null
#SBATCH        --output=/dev/null

# Extract chrom and pos of the target window
CHROM=$( sed -n $((${SLURM_ARRAY_TASK_ID}))p windows.txt | awk '{print $1}' )
POS=$( sed -n $((${SLURM_ARRAY_TASK_ID}))p windows.txt | awk '{print $2}' )
POS1=$((${POS}-10000))
POS2=$((${POS}+10000))

# Source required programs
. /local/env/envconda.sh
conda activate ~/TOOLS/angsd-0.921

# Output folder
output=/scratch/tcapblancq/DILS

# Reference genome
ref=~/Study-Mapping-WGS/genome/BST1.fa

# Create a file with information requires in the header of DILS
#touch ${output}/names.tmp
#for i in $(cat all_bam.list) 
#do 
#	name=${i/.final.bam/} 
#	echo `basename $name` >> ${output}/names.tmp 
#done
#paste -d "|" species.txt ${output}/names.tmp > ${output}/header.tmp
#cat ${output}/header.tmp | tr '\n' ' ' | sed -e 's/ /\t/g' > ${output}/header.txt
#rm ${output}/*.tmp

# Produce a geno file with genotypes for all individuals and positions (including monomorphic) for the genomic window 
angsd -b all_bam.list \
    -ref ${ref} \
    -anc ${ref} \
    -out ${output}/window${CHROM}:${POS1}-${POS2} \
    -r ${CHROM}:${POS1}-${POS2} \
    -nThreads 10 \
    -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
    -GL 1 \
    -doMaf 1 -doMajorMinor 1 \
    -doPost 1 -postCutoff 0.7 -doGeno 4 -doCounts 1 -geno_minDepth 3 -geno_maxDepth 200  

# Check the number of sites genotypes
nb_loci=$(zcat ${output}/window${CHROM}:${POS1}-${POS2}.geno.gz | wc -l)

if (("${nb_loci}" < "5000")) 
then 

	echo "Not enough sites on ${CHROM}:${POS1}-${POS2}" 
	
else

	# Format the file
  	zcat ${output}/window${CHROM}:${POS1}-${POS2}.geno.gz | cut -f3- > ${output}/window${CHROM}:${POS1}-${POS2}.tmp
  	cat ${output}/header.txt <(echo) ${output}/window${CHROM}:${POS1}-${POS2}.tmp > ${output}/geno_window${CHROM}:${POS1}-${POS2}.tmp

	# Some of Fred's magic !
	cat ${output}/geno_window${CHROM}:${POS1}-${POS2}.tmp | gawk '{gsub(/\t$/,"");}1' | gawk -F"\t" 'FNR==1{for (i=1;i<=NF;i++) {header[i] = $i;}nbInd=NF;next;}{for (i=1;i<=nbInd;i++) {seq[i,1] = seq[i,1] "" substr($i,1,1);seq[i,2] = seq[i,2] "" substr($i,2,1);}}END{for (i=1;i<=nbInd;i++) {for (j=1;j<=2;j++) print ">"header[i]"|allele"j"\n" seq[i,j]}}' > ${output}/window${CHROM}:${POS1}-${POS2}.dils
  	
	# Complete the header for each allele
  	sed -i "s/>/>window_${CHROM}:${POS}|/g" ${output}/window${CHROM}:${POS1}-${POS2}.dils

  	# Remove temp files
	rm ${output}/window${CHROM}:${POS1}-${POS2}.tmp
	rm ${output}/geno_window${CHROM}:${POS1}-${POS2}.tmp

fi

# Remove intermediate files
rm ${output}/window${CHROM}:${POS1}-${POS2}.geno.gz
rm ${output}/window${CHROM}:${POS1}-${POS2}.mafs.gz
rm ${output}/window${CHROM}:${POS1}-${POS2}.arg

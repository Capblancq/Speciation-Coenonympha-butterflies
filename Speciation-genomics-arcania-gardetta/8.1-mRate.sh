#!/bin/bash
#SBATCH -p genscale
#SBATCH --job-name=mRate_windows
#SBATCH --cpus-per-task=10
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --array=9001-9844%20
#SBATCH        --error=/dev/null
#SBATCH        --output=/dev/null

# Extract chrom and pos of the target window
CHROM=$( sed -n $((${SLURM_ARRAY_TASK_ID}))p windows.txt | awk '{print $1}' )
POS=$( sed -n $((${SLURM_ARRAY_TASK_ID}))p windows.txt | awk '{print $2}' )
POS1=$((${POS}-25000))
POS2=$((${POS}+25000))

# Input/Output folders
INPUT=/scratch/tcapblancq/Genotypes
mkdir /scratch/tcapblancq/mRate
OUTPUT=/scratch/tcapblancq/mRate

# Source ANGSD
. /local/env/envconda.sh
conda activate ~/TOOLS/angsd-0.921

# Estimates 2D-SFS from saf files
realSFS ${INPUT}/arcania_allsites.saf.idx ${INPUT}/gardetta_allsites.saf.idx -r ${CHROM}:${POS1}-${POS2} > ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}.sfs

# Source R
. /local/env/envr-3.5.1.sh

# Creates 2D-SFS with Fastsimcoal format
Rscript --vanilla 8.2-2D-SFS.R ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}.sfs ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}_jointMAFpop1_0.obs

# Add tabulation first line
sed -i 's/d0_0/\td0_0/g' ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}_jointMAFpop1_0.obs

# Add header
echo "1 observations" > ${OUTPUT}/header.txt
cat ${OUTPUT}/header.txt  ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}_jointMAFpop1_0.obs > ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}_jointMAFpop1_0.tmp
mv ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}_jointMAFpop1_0.tmp ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}_jointMAFpop1_0.obs

# Goes to working directory
cd ${OUTPUT}

# Creates directory for the scenario
mkdir window${SLURM_ARRAY_TASK_ID}

# Copies required files
cp SC_windows_m0.tpl SC_windows_m0.est window${SLURM_ARRAY_TASK_ID}/
cp ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}_jointMAFpop1_0.obs window${SLURM_ARRAY_TASK_ID}/SC_windows_m0_jointMAFpop1_0.obs
cp SC_windows_mfixed.tpl SC_windows_mfixed.est window${SLURM_ARRAY_TASK_ID}/
cp ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}_jointMAFpop1_0.obs window${SLURM_ARRAY_TASK_ID}/SC_windows_mfixed_jointMAFpop1_0.obs

# Enters the run directory
cd window${SLURM_ARRAY_TASK_ID}/

# Run fastsimcaol with varying me
~/TOOLS/fsc2705 -t SC_windows_m0.tpl -n1000000 -m -e SC_windows_m0.est -M -L40 -q -c10

# Run fastsimcaol with fixed me
~/TOOLS/fsc2705 -t SC_windows_mfixed.tpl -n1000000 -m -e SC_windows_mfixed.est -M -L40 -q -c10

# Estimate delta B following Laetsch et al. 2022
lnCL1=`cat ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}/SC_windows_m0/SC_windows_m0.bestlhoods | awk 'NR==2 {print $3}'`
lnCL2=`cat ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}/SC_windows_mfixed/SC_windows_mfixed.bestlhoods | awk 'NR==2 {print $3}'`
DeltaB=$(bc <<< "${lnCL1} - ${lnCL2}")

# Appends results to the results file
if test -f ${OUTPUT}/results_SC_m0_windows.txt
then
	echo "${CHROM} ${POS}" `cat ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}/SC_windows_m0/SC_windows_m0.bestlhoods | awk 'NR==2'` ${DeltaB} >> ${OUTPUT}/results_SC_m0_windows.txt
else
	echo "Chrom midPos" `cat ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}/SC_windows_m0/SC_windows_m0.bestlhoods | awk 'NR==1'` "DeltaB" > ${OUTPUT}/results_SC_m0_windows.txt
	echo ${CHROM} ${POS} `cat ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}/SC_windows_m0/SC_windows_m0.bestlhoods | awk 'NR==2'` ${DeltaB} >> ${OUTPUT}/results_SC_m0_windows.txt
fi

# Remove temporary files
rm -rf ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}/
rm ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}.sfs
rm ${OUTPUT}/window${SLURM_ARRAY_TASK_ID}_jointMAFpop1_0.obs

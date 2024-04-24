#!/bin/bash

##SBATCH -p genscale
#SBATCH --job-name=snpEff
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=30G
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null


## Source required programs
. /local/env/envconda.sh
conda activate ~/TOOLS/java 

## Output folder
output=/scratch/tcapblancq/Genotypes

## Building a new genome database for snpEff
#nano ~/TOOLS/snpEff/snpEff.config
## Adding 
#    # Coenonympha arcania
#    Carcania.genome : Coenonympha_arcania
#mkdir ~/TOOLS/snpEff/data/Carcania
#cp /groups/divalps/Annotation/Helixer/BST1/BST1_nucl_mito.fa ~/TOOLS/snpEff/data/Carcania/sequences.fa # genome
#cp /groups/divalps/Annotation/Helixer/BST1/CARC.gff3 ~/TOOLS/snpEff/data/Carcania/genes.gff # annotations
#cp /groups/divalps/Annotation/Helixer/BST1/CARC_cds.fa ~/TOOLS/snpEff/data/Carcania/cds.fa
#cp /groups/divalps/Annotation/Helixer/BST1/CARC_proteins.fa ~/TOOLS/snpEff/data/Carcania/protein.fa
#cp /groups/divalps/Annotation/Helixer/BST1/CARC_transcripts.fa ~/TOOLS/snpEff/data/Carcania/transcripts.fa

## Build the new dataset/reference for C. arcania
#java -jar ~/TOOLS/snpEff/snpEff.jar build -gff3 -v Carcania

# Using a vcf file with all the polymorphic position from ANGSD 
java -Xmx30G -jar ~/TOOLS/snpEff/snpEff.jar -c ~/TOOLS/snpEff/snpEff.config -v Carcania ${output}/All_polymorphicsites.vcf.gz > ${output}/All_polymorphicsites_annotated.vcf

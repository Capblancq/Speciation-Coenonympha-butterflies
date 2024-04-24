#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=30G
##SBATCH        --error=/dev/null
##SBATCH        --output=/dev/null

# Working Directory
mkdir /scratch/tcapblancq/Dxy
DIR=/scratch/tcapblancq/Dxy

# Allele frequency
#zcat /scratch/tcapblancq/Genotypes/arcania_allsites.mafs.gz | awk '{print $1,$2,$7}' > ${DIR}/arcania_mafs.txt 
#zcat /scratch/tcapblancq/Genotypes/gardetta_allsites.mafs.gz | awk '{print $1,$2,$7}' > ${DIR}/gardetta_mafs.txt

# Copy window file to the right directory 
#cp windows.txt ${DIR}/

# Use R to estimate Dxy for each windows of windows.txt
. /local/env/envr-3.5.1.sh

cd ${DIR}

R --vanilla << "EOF"

  # Data
  allfreqA <- read.table("arcania_mafs.txt", header=T)
  allfreqB <- read.table("gardetta_mafs.txt", header=T)   
  allfreq <- merge(allfreqA, allfreqB, by=c("chromo","position"))
  allfreq <- allfreq[order(allfreq$chromo, allfreq$position),]

  # Estimating Dxy from allele frequencies and for each site
  allfreq$dxy <- (allfreq$knownEM.x*(1-allfreq$knownEM.y))+(allfreq$knownEM.y*(1-allfreq$knownEM.x))
  write.table(allfreq[,c("chromo","position","dxy")], file="Dxy_persite.txt",quote=FALSE, row.names=FALSE, sep='\t')

  # Windows file
  win <- read.table("windows.txt", header = F, stringsAsFactors = F)

  # Dxy per window
  win$Dxy <- unlist(lapply(1:nrow(win), function(x) mean(allfreq$dxy[which(as.character(allfreq$chromo)==as.character(win[x,1]) & allfreq$position>=as.numeric(win[x,2]-25000) & allfreq$position<as.numeric(win[x,2]+25000))])))
  write.table(win, file="Dxy_perwindow.txt",quote=FALSE, row.names=FALSE, sep='\t')

EOF

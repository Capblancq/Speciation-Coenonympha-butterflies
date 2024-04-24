#!/usr/bin/env Rscript

# Arguments
args = commandArgs(trailingOnly=TRUE)

# 2D SFS arcania (pop0) / gardetta (pop1)
sfs <- scan(paste(args[1], sep=""), quiet=T)
N1 <- 19
N2 <- 16
tab <- matrix(data = sfs, ncol = (N1*2+1), nrow = (N2*2+1))
colnames(tab) <- paste(rep("d0", (N1*2+1)), seq(0,(N1*2)), sep = "_")
row.names(tab) <- paste(rep("d1", (N2*2+1)), seq(0,(N2*2)), sep = "_")

# Function to transform dadi like SFS into 2D-SFS for Fastsimcoal2
derived2maf=function(derived_sfs2d){
	 n1=nrow(derived_sfs2d)
	 n2=ncol(derived_sfs2d)
	 maf_2dsfs=matrix(0,nrow=n1,ncol=n2)
	 colnames(maf_2dsfs)=colnames(derived_sfs2d)
	 rownames(maf_2dsfs)=rownames(derived_sfs2d)
	 threshold_freq=0.5*(n1+n2-2)
	 for(i in 0:(n1-1)){
	   for(j in 0:(n2-1)){
	     if(i+j < threshold_freq){
	       maf_2dsfs[i+1,j+1]=maf_2dsfs[i+1,j+1]+derived_sfs2d[i+1,j+1]
	     }
	     else if(i+j == threshold_freq){
	       maf_2dsfs[i+1,j+1]=maf_2dsfs[i+1,j+1]+(0.5*derived_sfs2d[i+1,j+1]+0.5*derived_sfs2d[n1-i,n2-j])
	     }
	     else{
	       maf_2dsfs[n1-i,n2-j]=maf_2dsfs[n1-i,n2-j]+derived_sfs2d[i+1,j+1]
	     }
	   }
	 }
	 maf_2dsfs
}

# Convert SFS 
tab2 <- derived2maf(tab)

# Write output with the right name for Fastsimcoal2
write.table(tab2, args[2], row.names = T, col.names = T, sep = "\t", quote = F)

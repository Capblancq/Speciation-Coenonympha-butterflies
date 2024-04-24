
cd /scratch/tcapblancq/Fst

. /local/env/envr-3.5.1.sh

R --vanilla

	# 2D SFS arcania (pop0) / gardetta (pop1)
	sfs <- scan(paste("arcania.gardetta.sfs", sep=""), quiet=T)
	N1 <- 19
	N2 <- 16
	tab <- matrix(data = sfs, ncol = (N1*2+1), nrow = (N2*2+1))
	colnames(tab) <- paste(rep("d0", (N1*2+1)), seq(0,(N1*2)), sep = "_")
	row.names(tab) <- paste(rep("d1", (N2*2+1)), seq(0,(N2*2)), sep = "_")
	write.table(tab, "Coenonympha_jointDAFpop1_0.obs", row.names = T, col.names = T, sep = "\t", quote = F)

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

	tab2 <- derived2maf(tab)
	write.table(tab2, "Coenonympha_jointMAFpop1_0.obs", row.names = T, col.names = T, sep = "\t", quote = F)

EOF

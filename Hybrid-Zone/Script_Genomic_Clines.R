##############################################################################################
##############################################################################################

# Script genetic clines analyses (Coenonypmha) 
# Capblancq Thibaut
# Date : 31-10-2018

##############################################################################################
##############################################################################################

library(introgress)
library(raster)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(cowplot)
library(HIest)
library(hzar)
library(vcfR)
library(adegenet)
library(foreach)
library(doParallel)


########################
#### Data filtering ####

# Information about samples
info<-read.table("infos_ind_VARS", header=T, sep = ",")

# Genetic data
GenData <- read.table("Coenonympha_1SNPrandom.gen", header = T)
row.names(GenData)<- GenData$Ind

################
#### HINDEX ####

loci.data<-as.data.frame(cbind(colnames(GenData)[3:ncol(GenData)],rep("C",length(colnames(GenData)[3:ncol(GenData)]))))
colnames(loci.data)<-c("Locus","type")

## Formating data
data_VARS<-prepare.data(admix.gen = t(GenData[which(info$Zone=="VARS"),]), loci.data = loci.data, parental2 = t(GenData[which(info$Zone=="GARDETTA"),-c(1:2)]), parental1 = t(GenData[which(info$Zone=="MACROMMA"),-c(1:2)])) 

## Heterozygosity
het_VARS<-calc.intersp.het(data_VARS)

## HINDEX
HINDEX_VARS<-est.h(introgress.data = data_VARS, loci.data = loci.data, ind.touse = NULL, fixed = FALSE)
hist(HINDEX_VARS$h, breaks = seq(0,1,0.1), col="darkgrey", xlim=c(0,1))
row.names(HINDEX_VARS)<-GenData[which(info$Zone=="VARS"),"Ind"]

## Heterozygotie ~ HINDEX
pdf("./Het_HINDEX.pdf", width = 4.5, height = 4.5)
plot(HINDEX_VARS$h, het_VARS, ylim = c(0,1), xlim = c(0,1), xlab = "HINDEX", ylab = "Heterozygosity", pch = 21, col= "black", bg="gray49")
dev.off()

## SNPs divergent in parental populations
tab<-cbind(as.data.frame(data_VARS$Parental1.allele.freq), as.data.frame(data_VARS$Parental2.allele.freq))
snps_fixed<-c(row.names(tab[which(tab[,1]>=0.9 & tab[,4]>=0.9),]), row.names(tab[which(tab[,2]>=0.9 & tab[,3]>=0.9),]))

## And with only one SNPs by ddRAD tag
snps_fixed <- as.vector(aggregate(snps_fixed, by = list(unlist(lapply(strsplit(snps_fixed, split = "_"), function(x) x[1]))), FUN = head, 1)[,2])

## Geographic clines 
dist_Vars<-pointDistance(c(997027.1, 6383629.52),info[info$Zone=="VARS",4:5], lonlat=FALSE)
pdf("./HINDEX.pdf", width = 5, height = 5)
plot(HINDEX_VARS$h~dist_Vars, ylim=c(0,1), cex=1.5, pch = 21, col= "black", bg="gray49", xlab= "Distance along hybrid zone (m)", ylab="Hybrid index")
dev.off()


########################
#### Genetic clines ####

### Frequences alleliques par site 
countpop<-NULL
nbindpop<-NULL
freq<-NULL
Popfreq<-list()
for (i in info$Pop[info$Zone=="VARS"])
{countpop<-apply(data_VARS$Count.matrix[snps_fixed,which(info$Pop[info$Zone=="VARS"]==i)], 1,sum, na.rm=T)
nbindpop<-ncol(data_VARS$Count.matrix[snps_fixed,which(info$Pop[info$Zone=="VARS"]==i)])-apply(apply(data_VARS$Count.matrix[snps_fixed,which(info$Pop[info$Zone=="VARS"]==i)], 1, is.na), 2, sum)
freq<-data.frame(A1=countpop/(nbindpop*2), A2=(1-countpop/(nbindpop*2)), NbInd=nbindpop)
Popfreq[[i]]<-freq
}

dist_Vars_pop <- aggregate(dist_Vars, by = list(info$Pop[info$Zone=="VARS"]), mean)
Popfreq <- Popfreq[as.character(dist_Vars_pop$Group.1[order(dist_Vars_pop$x,decreasing = FALSE)])]
dist_Vars_pop <- dist_Vars_pop[order(dist_Vars_pop$x,decreasing = FALSE),]

PopFreqA1<-do.call(data.frame, lapply(Popfreq, function(Popfreq) Popfreq[,1]))
row.names(PopFreqA1)<-snps_fixed

PopFreqA2<-do.call(data.frame, lapply(Popfreq, function(Popfreq) Popfreq[,2]))
row.names(PopFreqA2)<-snps_fixed

NbInd<-do.call(data.frame, lapply(Popfreq, function(Popfreq) Popfreq[,3]))
row.names(NbInd)<-snps_fixed

freqglob<-apply(PopFreqA2, 2, mean, na.rm=T)
nbind<-apply(NbInd, 2, max, na.rm=T)

hindexglob<-aggregate(HINDEX_VARS$h, by = list(info$Pop[info$Zone=="VARS"]), mean)
row.names(hindexglob) <- hindexglob$Group.1
hindexglob <- hindexglob[c("TOU", "SPA", "SPA_2","SPA_3", "MEL_INF", "MEL_SUP", "SLC", "VAR_1", "VAR_2", "VAR_3", "VAR_N"),"x"]

### Mise en forme des donnees pour la frequence moyenne
c.snps_VARS<-lapply(1:nrow(PopFreqA2), function(x) hzar.doMolecularData1DPops(distance = dist_Vars_pop$x, pObs = as.numeric(PopFreqA2[x,]), nEff=as.numeric(nbind)))
c.snps_VARS[[length(c.snps_VARS)+1]]<-hzar.doMolecularData1DPops(distance = dist_Vars_pop$x, pObs = as.numeric(freqglob), nEff=as.numeric(nbind))
c.snps_VARS[[length(c.snps_VARS)+1]]<-hzar.doMolecularData1DPops(distance = dist_Vars_pop$x, pObs = as.numeric(hindexglob), nEff=as.numeric(nbind))
names(c.snps_VARS)<-c(row.names(PopFreqA2), "freqglobal", "hindex")

### Fitting clines
hzar_function <- function(TABLE)
{
  comb1<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "fixed", tails="none", direction = NULL)
  firstfit_comb1<-hzar::hzar.first.fitRequest.old.ML(model=comb1, obsData=TABLE)
  fit_comb1<-hzar::hzar.chain.doSeq(firstfit_comb1, 3, collapse=T)
  
  comb2<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "fixed", tails="right", direction = NULL)
  firstfit_comb2<-hzar::hzar.first.fitRequest.old.ML(model=comb2, obsData=TABLE)
  fit_comb2<-hzar::hzar.chain.doSeq(firstfit_comb2, 3, collapse=T)
  
  comb3<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "fixed", tails="left", direction = NULL)
  firstfit_comb3<-hzar::hzar.first.fitRequest.old.ML(model=comb3, obsData=TABLE)
  fit_comb3<-hzar::hzar.chain.doSeq(firstfit_comb3, 3, collapse=T)
  
  comb4<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "fixed", tails="mirror", direction = NULL)
  firstfit_comb4<-hzar::hzar.first.fitRequest.old.ML(model=comb4, obsData=TABLE)
  fit_comb4<-hzar::hzar.chain.doSeq(firstfit_comb4, 3, collapse=T)
  
  comb5<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "none", tails="both", direction = NULL)
  firstfit_comb5<-hzar::hzar.first.fitRequest.old.ML(model=comb5, obsData=TABLE)
  fit_comb5<-hzar::hzar.chain.doSeq(firstfit_comb5, 3, collapse=T)
  
  comb6<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "none", tails="none", direction = NULL)
  firstfit_comb6<-hzar::hzar.first.fitRequest.old.ML(model=comb6, obsData=TABLE)
  fit_comb6<-hzar::hzar.chain.doSeq(firstfit_comb6, 3, collapse=T)
  
  comb7<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "none", tails="right", direction = NULL)
  firstfit_comb7<-hzar::hzar.first.fitRequest.old.ML(model=comb7, obsData=TABLE)
  fit_comb7<-hzar::hzar.chain.doSeq(firstfit_comb7, 3, collapse=T)
  
  comb8<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "none", tails="left", direction = NULL)
  firstfit_comb8<-hzar::hzar.first.fitRequest.old.ML(model=comb8, obsData=TABLE)
  fit_comb8<-hzar::hzar.chain.doSeq(firstfit_comb8, 3, collapse=T)
  
  comb9<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "none", tails="mirror", direction = NULL)
  firstfit_comb9<-hzar::hzar.first.fitRequest.old.ML(model=comb9, obsData=TABLE)
  fit_comb9<-hzar::hzar.chain.doSeq(firstfit_comb9, 3, collapse=T)
  
  comb10<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "none", tails="both", direction = NULL)
  firstfit_comb10<-hzar::hzar.first.fitRequest.old.ML(model=comb10, obsData=TABLE)
  fit_comb10<-hzar::hzar.chain.doSeq(firstfit_comb10, 3, collapse=T)
  
  comb11<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "free", tails="none", direction = NULL)
  firstfit_comb11<-hzar::hzar.first.fitRequest.old.ML(model=comb11, obsData=TABLE)
  fit_comb11<-hzar::hzar.chain.doSeq(firstfit_comb11, 3, collapse=T)
  
  comb12<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "free", tails="right", direction = NULL)
  firstfit_comb12<-hzar::hzar.first.fitRequest.old.ML(model=comb12, obsData=TABLE)
  fit_comb12<-hzar::hzar.chain.doSeq(firstfit_comb12, 3, collapse=T)
  
  comb13<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "free", tails="left", direction = NULL)
  firstfit_comb13<-hzar::hzar.first.fitRequest.old.ML(model=comb13, obsData=TABLE)
  fit_comb13<-hzar::hzar.chain.doSeq(firstfit_comb13, 3, collapse=T)
  
  comb14<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "free", tails="mirror", direction = NULL)
  firstfit_comb14<-hzar::hzar.first.fitRequest.old.ML(model=comb14, obsData=TABLE)
  fit_comb14<-hzar::hzar.chain.doSeq(firstfit_comb14, 3, collapse=T)
  
  comb15<-hzar::hzar.makeCline1DFreq(data=TABLE, scaling = "free", tails="both", direction = NULL)
  firstfit_comb15<-hzar::hzar.first.fitRequest.old.ML(model=comb15, obsData=TABLE)
  fit_comb15<-hzar::hzar.chain.doSeq(firstfit_comb15, 3, collapse=T)
  
  reslist <- list(fit_comb1, fit_comb2, fit_comb3, fit_comb4, fit_comb5, fit_comb6, fit_comb7, fit_comb8, fit_comb9, fit_comb10, fit_comb11, fit_comb12, fit_comb13, fit_comb14, fit_comb15)
  
  return(reslist)
}

#setup parallel backend to use many processors
registerDoParallel(cores = 12)
list_clines_snp_VARS <- list()
list_clines_snp_VARS <- foreach(i=c.snps_VARS) %dopar% {hzar_function(i)}

group_snp_VARS<-lapply(list_clines_snp_VARS, hzar.make.obsDataGroup)


### Model selection with AIC criterium
c.snp_VARS<-list()
modelid<-NULL
for (i in 1:length(c.snps_VARS))
{aic<-hzar.AIC.hzar.obsDataGroup(group_snp_VARS[[i]])
c.snp_VARS[[i]]<-hzar.make.obsDataGroup(list_clines_snp_VARS[[i]][[which.min(as.vector(aic$AIC))]])
modelid<-c(modelid,which.min(as.vector(aic$AIC)))
}

### Center and slope of the clines
c<-NULL
w<-NULL
width<-NULL
center<-NULL
mm<-NULL
pminmax<-NULL
par<-NULL
params_VARS<-NULL
for (i in 1:length(c.snp_VARS))
{
  c<-density(list_clines_snp_VARS[[i]][[which.min(as.vector(aic$AIC))]][[1]][[5]][,1])
  center<-c$x[which.max(c$y)]
  w<-density(list_clines_snp_VARS[[i]][[which.min(as.vector(aic$AIC))]][[1]][[5]][,2])
  width<-w$x[which.max(w$y)]
  toto<-c.snp_VARS[[i]]
  mm<-hzar.get.ML.cline(toto$data.groups[[1]])
  pminmax<-c(mm$param.all$pMin, mm$param.all$pMax)
  par<-c(center, width, pminmax)
  params_VARS<-rbind(params_VARS, par)
}
row.names(params_VARS)<-names(c.snps_VARS)
colnames(params_VARS)<-c("center", "width", "Pmin", "Pmax")
params_VARS<-as.data.frame(params_VARS)
params_VARS$slope<-(params_VARS$Pmax-params_VARS$Pmin)/params_VARS$width

deltaP<-hist(params_VARS$Pmax-params_VARS$Pmin, breaks=seq(0,1,0.05))[2][[1]]

pdf("DeltaFreq.pdf", width = 5, height = 5)
barplot(deltaP, col="gray49", xlab = "Pmax-Pmin", ylab= "Count")
dev.off()

pdf("CenterSlope.pdf", width = 5, height = 5)
plot(params_VARS$slope~params_VARS$center, xlim=c(0, 8000), col= "black", bg="gray49", pch=21, xlab = "Center of the cline", ylab = "Slope of the cline")
points(params_VARS$center[length(params_VARS$center)], params_VARS$slope[length(params_VARS$center)], col= "black", bg="red", pch=24, cex=1.5)
dev.off()


## Visualization

## Little change in the plot.cline function of hzar
hzar.plot.cline2 <- function(cline, add = FALSE, ylim = FALSE, ...) 
{
  x = NULL
  if (inherits(cline, "hzar.cline")) 
    curve(cline$clineFunc(x), add = add, ...)
  if (inherits(cline, c("hzar.dataGroup", "hzar.fitRequest"))) {
    dataGroup <- hzar.fit2DataGroup(cline)
    #hzar.plot.obsData(dataGroup, add = add, ylim = ylim, ...)
    hzar.plot.cline(hzar.get.ML.cline(dataGroup), add = TRUE, ...)
  }
  if (inherits(cline, c("hzar.obsDataGroup"))) {
    if (identical(ylim, FALSE)) {
      yS <- sapply(cline$data.groups, function(x) {
        ylim <- c(0, 1)
        try(ylim <- x$obsData$ylim)
        as.numeric(ylim)
      })
      ylim <- c(min(yS), max(yS))
    }
    hzar.plot.obsData(cline, add = add, ylim = ylim, cex=0, ...)
    lapply(cline$data.groups, function(dataGroup) hzar.plot.cline(hzar.get.ML.cline(dataGroup), add = TRUE, ...))
  }
}


pdf("genetic_clines_Vars.pdf", width = 5, height = 5)
hzar.plot.cline2(c.snp_VARS[[length(c.snp_VARS)]], col="white", add=FALSE)
for (i in 1:(length(c.snp_VARS)-2))
{hzar.plot.cline2(c.snp_VARS[[i]], add = TRUE, col="gray34", lwd=1.5)}
#hzar.plot.cline2(c.snp_VARS[[length(c.snp_VARS)-1]], col="red", add=TRUE, lwd=4, lty=2)
hzar.plot.cline2(c.snp_VARS[[length(c.snp_VARS)]], col="red", add=TRUE, lwd=4, lty=2)
dev.off()

## Table cline parameters
write.table(params_VARS, "./Cline_parameters.csv")

## SNP with high slope parameters
snps_fixed[which(params_VARS$slope>0.0005)]
# "102796_29" "109431_31" "112938_28" "1159_23"   "1160_43"   "143705_7" "144664_17" "1630_29"   "1650_49"   "1825_31"   "1941_71"   "2092_16"  "2129_94"   "2467_94"   "3259_23"   "3620_61"   "3650_11"   "4194_93"  "4203_12"   "4520_52"   "5199_70"   "7736_71"   "851_56"



##############################################################################################
##############################################################################################

# Script pour calculer les indices d'hybridation Hindex (Coenonypmha) 
# Capbalancq Thibaut
# Date : 15-09-2015

##############################################################################################
##############################################################################################

setwd("/Users/capblancq/Documents/Projet Papillons/These Hybridation Coenonympha/ddRADseq/Analyses ddRADseq/HINDEX/Indice_Hybridation/")

##### Calcul du HINDEX, de l'heterozygotie et mise en forme graphique (Buerkle 2005)

library(introgress)
library(raster)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(cowplot)

#### Chargement des donnees totales
geno<-read.table("Input_Introgress_390ind_r0.4p4maf0.01.txt", header=T, check.names = FALSE) 
Loca<-NULL
for (i in 1:nrow(geno)){Loca<-c(Loca,strsplit(as.character(geno[i,"Ind"]),"_")[[1]][1])}
Loc<-NULL
for (i in 1:length(Loca)){Loc<-c(Loc,strsplit(as.character(Loca[i]),"-")[[1]][1])}

## Pour differentier les darwiniana des arcania a chaillol et orciere et les gardetta des arcania à Vaujany
#test<-cbind(Loc[c(13:15, 304:323, 184:205, 263:286, 344:346)],c(rep("LNO_dar",23), rep("VJ_gar",22),rep("CHL_dar",24),rep("ORC_dar",3)))
#Loc[c(13:15, 304:323, 184:205, 263:286, 344:346)]<-c(rep("LNO_dar",23),rep("VJ_gar",22),rep("CHL_dar",24),rep("ORC_dar",3))

loci.data<-as.data.frame(cbind(colnames(geno)[3:ncol(geno)],rep("C",length(colnames(geno)[3:ncol(geno)])))) ## Informations sur la nature des loci = co-dominants dans notre cas
colnames(loci.data)<-c("Locus","type") ## Probleme dans la fonction est.h il demande la colonne type donc on est obligé de l'appeler comme ça, peut-être que ça serait mieux des le debut...

##############################################
##### Pour les individus de Vars/St-Paul #####

## Mise en forme des donnees
data_vars<-prepare.data(admix.gen = t(geno[which(Loc=="SLC"|Loc=="MEL"),]), loci.data = as.data.frame(cbind(colnames(geno)[3:ncol(geno)],rep("C",length(colnames(geno)[3:ncol(geno)])))), parental2 = t(geno[which(Loc=="LAU"|Loc=="SES"|Loc=="MOO"|Loc=="AIL"|Loc=="SEL"),3:ncol(geno)]), parental1 = t(geno[which(Loc=="LAR"|Loc=="BOR"|Loc=="FOA"|Loc=="LOM"|Loc=="SEY"),3:ncol(geno)])) ## Pour preparer les donnees avec uniquement les individus de SLC, VAR, SPA et MEL comme hybrides potentiels et les gardetta et darwiniana proches comme parentaux

## Calcul de l'heterozygotie chez les individus/pops de la zone hybride
het_vars<-calc.intersp.het(data_vars) ## Calcul de l'heterozygotie pour les individus de la pop hybride

## Calcul du HINDEX chez les individus de la zone hybride
HINDEX_vars<-est.h(introgress.data = data_vars, loci.data = loci.data, ind.touse = NULL, fixed = FALSE, p1.allele=NULL, p2.allele=NULL)
hist(HINDEX_vars$h, breaks = seq(0,1,0.1), col="darkgrey", xlim=c(0,1))

row.names(HINDEX_vars)<-geno[which(Loca=="SLC"|Loca=="MEL"),"Ind"]

## HINDEX en fonction de l'heterozygotie
triangle_plot_bis(HINDEX_vars, het_vars, out.file="tri_plot_Vars.pdf") ## Permet d'estimer graphiquement si les hybrides sont des F1, backcross... Ou si la pop hybride est plutôt stable avecune heterozygotie pas vraiment superieure a celle des parents

## Clines génomiques
cline_out_vars<- genomic.clines(data_vars, HINDEX_vars, loci.data, sig.test=T, method="parametric", n.reps=100, classification=T)
clines.plot(cline_out_vars, out.file="cline_plot_Vars.pdf")

## Marker ancestry across individuals

mk.image(data_vars, loci.data, marker.order=NULL, HINDEX_vars, ind.touse=NULL, loci.touse=NULL, pdf=FALSE)


##############################################
##### Idem pour la zone de contact Ecrins ####

data_ecrins<-prepare.data(admix.gen = t(geno[which(Loc=="AIL"|Loc=="NAV"|Loc=="GIO"|Loc=="CHL_dar"|Loc=="ORC_dar"|Loc=="BY"|Loc=="VALS"|Loc=="SEL"|Loc=="DOR"|Loc=="NDS_dargar"),]), loci.data = as.data.frame(cbind(colnames(geno)[3:ncol(geno)],rep("C",length(colnames(geno)[3:ncol(geno)])))), parental1 = t(geno[which(Loc=="LAU"|Loc=="SES"|Loc=="MOO"|Loc=="VJ_gar"),3:ncol(geno)]), parental2 = t(geno[which(Loc=="LAR"|Loc=="BOR"|Loc=="FOA"|Loc=="LOM"|Loc=="SEY"),3:ncol(geno)])) ## Pour preparer les donnees avec uniquement les individus de SLC, VAR, SPA et MEL comme hybrides potentiels et les gardetta et darwiniana proches comme parentaux

het_ecrins<-calc.intersp.het(data_ecrins) 

HINDEX_ecrins<-est.h(introgress.data = data_ecrins, loci.data = loci.data, ind.touse = NULL, fixed = FALSE, p1.allele=NULL, p2.allele=NULL)
hist(HINDEX_ecrins$h)

row.names(HINDEX_ecrins)<-geno[which(Loc=="AIL"|Loc=="NAV"|Loc=="GIO"|Loc=="CHL_dar"|Loc=="ORC_dar"|Loc=="BY"|Loc=="VALS"|Loc=="SEL"|Loc=="DOR"),"Ind"]

triangle.plot(HINDEX_ecrins, het_ecrins, out.file="tri_plot_Ecrins.pdf") ## Permet d'estimer graphiquement si les hybrides sont des F1, backcross... Ou si la pop hybride est plutôt stable avecune heterozygotie pas vraiment superieure a celle des parents

cline_out_ecrins<- genomic.clines(data_ecrins, HINDEX_ecrins, loci.data, sig.test=T, method="parametric", n.reps=100, classification=T)
clines.plot(cline_out_ecrins, out.file="cline_plot_Ecrins.pdf")

mk.image(data_ecrins, loci.data, marker.order=NULL, HINDEX_ecrins, ind.touse=NULL, loci.touse=NULL, pdf=FALSE)

##############################################
##### Idem les individus de Bellwald/Moosalps #####

## Mise en forme des donnees
data_bel<-prepare.data(admix.gen = t(geno[which(Loc=="BEL"|Loc=="MOO"),]), loci.data = as.data.frame(cbind(colnames(geno)[3:ncol(geno)],rep("C",length(colnames(geno)[3:ncol(geno)])))), parental1 = t(geno[which(Loc=="LAU"|Loc=="SES"|Loc=="ALB"|Loc=="MOO"),3:ncol(geno)]), parental2 = t(geno[which(Loc=="ACQ"|Loc=="BSG"|Loc=="FTN"|Loc=="SIMD"),3:ncol(geno)])) ## Pour preparer les donnees avec uniquement les individus de SLC, VAR, SPA et MEL comme hybrides potentiels et les gardetta et darwiniana proches comme parentaux

## Calcul de l'heterozygotie chez les individus/pops de la zone hybride
het_bel<-calc.intersp.het(data_bel) ## Calcul de l'heterozygotie pour les individus de la pop hybride

## Calcul du HINDEX chez les individus de la zone hybride
HINDEX_bel<-est.h(introgress.data = data_bel, loci.data = loci.data, ind.touse = NULL, fixed = FALSE, p1.allele=NULL, p2.allele=NULL)
hist(HINDEX_bel$h, breaks = seq(0,1,0.1), col="darkgrey", xlim=c(0,1))

row.names(HINDEX_bel)<-geno[which(Loc=="BEL"|Loc=="MOO"),"Ind"]

## HINDEX en fonction de l'heterozygotie
triangle_plot_bis(HINDEX_bel, het_bel, out.file="tri_plot_Bel.pdf") ## Permet d'estimer graphiquement si les hybrides sont des F1, backcross... Ou si la pop hybride est plutôt stable avecune heterozygotie pas vraiment superieure a celle des parents

## Clines génomiques
cline_out_bel<- genomic.clines(data_bel, HINDEX_bel, loci.data, sig.test=T, method="parametric", n.reps=100, classification=T)
clines.plot(cline_out_bel, out.file="cline_plot_Bel.pdf")

#######################################################
##### Idem pour la zone de contact arc/gar Vaujany ####

data_vaujany<-prepare.data(admix.gen = t(geno[which(Loc=="VJ"),]), loci.data = as.data.frame(cbind(colnames(geno)[3:ncol(geno)],rep("C",length(colnames(geno)[3:ncol(geno)])))), parental1 = t(geno[which(Loc=="LAU"|Loc=="SES"|Loc=="AIL"|Loc=="SEL"),3:ncol(geno)]), parental2 = t(geno[which(Loc=="MJA"|Loc=="BAR"|Loc=="COL"|Loc=="VIG1"|Loc=="SDA"),3:ncol(geno)]))
het_vaujany<-calc.intersp.het(data_vaujany) 

HINDEX_vaujany<-est.h(introgress.data = data_vaujany, loci.data = loci.data, ind.touse = NULL, fixed = FALSE, p1.allele=NULL, p2.allele=NULL)
hist(HINDEX_vaujany$h, col="darkgrey")

triangle_plot_bis(HINDEX_vaujany, het_vaujany, out.file="tri_plot_Vaujany.pdf") ## Permet d'estimer graphiquement si les hybrides sont des F1, backcross... Ou si la pop hybride est plutôt stable avecune heterozygotie pas vraiment superieure a celle des parents

cline_out_vaujany<- genomic.clines(data_vaujany, HINDEX_vaujany, loci.data, sig.test=T, method="parametric", n.reps=100, classification=T)
clines.plot(cline_out_vaujany, out.file="cline_plot_Vaujany.pdf")

mk.image(data_vaujany, loci.data, marker.order=NULL, HINDEX_vaujany, ind.touse=NULL, loci.touse=NULL, pdf=FALSE)

#####################################################
##### Idem pour la zone de contact arc/gar Aoste ####

data_aoste<-prepare.data(admix.gen = t(geno[which(Loc=="AOS"),]), loci.data = as.data.frame(cbind(colnames(geno)[3:ncol(geno)],rep("C",length(colnames(geno)[3:ncol(geno)])))), parental1 = t(geno[which(Loc=="LAU"|Loc=="SES"|Loc=="MOO"|Loc=="AIL"|Loc=="SEL"|Loc=="OBE"),3:ncol(geno)]), parental2 = t(geno[which(Loc=="MJA"|Loc=="BAR"|Loc=="COL"|Loc=="VIG1"|Loc=="SDA"),3:ncol(geno)])) ## Pour preparer les donnees avec uniquement les individus de SLC, VAR, SPA et MEL comme hybrides potentiels et les gardetta et darwiniana proches comme parentaux

het_aoste<-calc.intersp.het(data_aoste) 

HINDEX_aoste<-est.h(introgress.data = data_aoste, loci.data = loci.data, ind.touse = NULL, fixed = FALSE, p1.allele=NULL, p2.allele=NULL)
hist(HINDEX_aoste$h, breaks = seq(0,1,0.1), col="darkgrey")
row.names(HINDEX_aoste)<-geno[which(Loc=="AOS"),"Ind"]

triangle_plot_bis(HINDEX_aoste, het_aoste, out.file="tri_plot_Aoste.pdf") ## Permet d'estimer graphiquement si les hybrides sont des F1, backcross... Ou si la pop hybride est plutôt stable avecune heterozygotie pas vraiment superieure a celle des parents

cline_out_aoste<- genomic.clines(data_aoste, HINDEX_aoste, loci.data, sig.test=T, method="parametric", n.reps=100, classification=T)
clines.plot(cline_out_aoste, out.file="cline_plot_Aoste.pdf")

########################################################
##### Idem pour la zone de contact arc/dar Chaillol ####

data_chaillol<-prepare.data(admix.gen = t(geno[which(Loc=="CHL"),]), loci.data = as.data.frame(cbind(colnames(geno)[3:ncol(geno)],rep("C",length(colnames(geno)[3:ncol(geno)])))), parental2 = t(geno[which(Loc=="LAR"|Loc=="BOR"|Loc=="FOA"|Loc=="LOM"|Loc=="SEY"),3:ncol(geno)]), parental1 = t(geno[which(Loc=="MJA"|Loc=="BAR"|Loc=="COL"|Loc=="VIG1"|Loc=="SDA"),3:ncol(geno)])) ## Pour preparer les donnees avec uniquement les individus de SLC, VAR, SPA et MEL comme hybrides potentiels et les gardetta et darwiniana proches comme parentaux

het_chaillol<-calc.intersp.het(data_chaillol) 

HINDEX_chaillol<-est.h(introgress.data = data_chaillol, loci.data = loci.data, ind.touse = NULL, fixed = FALSE, p1.allele=NULL, p2.allele=NULL)
hist(HINDEX_chaillol$h, col="darkgrey")

triangle_plot_bis(HINDEX_chaillol, het_chaillol, out.file="tri_plot_Chaillol.pdf") ## Permet d'estimer graphiquement si les hybrides sont des F1, backcross... Ou si la pop hybride est plutôt stable avecune heterozygotie pas vraiment superieure a celle des parents

cline_out_chaillol<- genomic.clines(data_chaillol, HINDEX_chaillol, loci.data, sig.test=T, method="parametric", n.reps=100, classification=T)
clines.plot(cline_out_chaillol, out.file="cline_plot_Chaillol.pdf")

########################################################
##### Idem pour la zone de contact arc/dar Locarno ####

data_locarno<-prepare.data(admix.gen = t(geno[which(Loc=="LNO"|Loc=="BAR"),]), loci.data = as.data.frame(cbind(colnames(geno)[3:ncol(geno)],rep("C",length(colnames(geno)[3:ncol(geno)])))), parental2 = t(geno[which(Loc=="ACQ"|Loc=="BSG"|Loc=="FTN"|Loc=="SIMD"),3:ncol(geno)]), parental1 = t(geno[which(Loc=="MJA"|Loc=="BAR"|Loc=="COL"|Loc=="VIG1"|Loc=="SDA"),3:ncol(geno)])) ## Pour preparer les donnees avec uniquement les individus de SLC, VAR, SPA et MEL comme hybrides potentiels et les gardetta et darwiniana proches comme parentaux

het_locarno<-calc.intersp.het(data_locarno) 

HINDEX_locarno<-est.h(introgress.data = data_locarno, loci.data = loci.data, ind.touse = NULL, fixed = FALSE, p1.allele=NULL, p2.allele=NULL)
hist(HINDEX_locarno$h, breaks = seq(0,1,0.1), col="darkgrey")

triangle_plot_bis(HINDEX_locarno, het_locarno, out.file="tri_plot_locarno.pdf") ## Permet d'estimer graphiquement si les hybrides sont des F1, backcross... Ou si la pop hybride est plutôt stable avecune heterozygotie pas vraiment superieure a celle des parents

cline_out_locarno<- genomic.clines(data_locarno, HINDEX_locarno, loci.data, sig.test=T, method="parametric", n.reps=100, classification=T)
clines.plot(cline_out_locarno, out.file="cline_plot_locarno.pdf")


#####################################

#### Figure HINDEX et hétérozygotie sur le même graphe #####

## Plot classique
pdf("HINDEX.pdf")
par(mfrow=c(3,2), mar= c(5, 4, 4, 4) + 0.3)
hist(HINDEX_vars$h, freq=TRUE, breaks = seq(0,1,0.1), col="grey88", xlim=c(0,1), xlab="HINDEX", las=1, ylab="Counts", main="vars")
#d <- density(HINDEX_vars$h)
#plot(d, add.plot=TRUE, xlim=c(0,1), xlab="HINDEX", main=NULL)
par(new=TRUE)
plot(HINDEX_vars$h, het_vars, axes=FALSE, xlab = "", ylab = "", bg="red", xlim=c(0,1), ylim=c(0,1), pch=21, cex=1.5, col="black")
axis(side=4, at = seq(0,1, 0.2), tick=TRUE, las=1, col="red", col.axis="red")
mtext("Heterozygosity", side=4, line=3, col="red")
dev.off()

## ggplot
table_vars<-data.frame(HINDEX=HINDEX_vars$h,Heterozygosity=het_vars)
table_bel<-data.frame(HINDEX=HINDEX_bel$h,Heterozygosity=het_bel)
table_aoste<-data.frame(HINDEX=HINDEX_aoste$h,Heterozygosity=het_aoste)
table_vaujany<-data.frame(HINDEX=HINDEX_vaujany$h,Heterozygosity=het_vaujany)
table_locarno<-data.frame(HINDEX=HINDEX_locarno$h,Heterozygosity=het_locarno)
table_chaillol<-data.frame(HINDEX=HINDEX_chaillol$h,Heterozygosity=het_chaillol)

points_aoste<-ggplot(table_aoste, aes(HINDEX, Heterozygosity)) + geom_point(size=6, position="jitter") + xlim(0,1) + ylim(0,0.5) + background_grid(major = 'y', minor = "none") + panel_border() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())
hist_aoste<-ggplot(table_aoste, aes(HINDEX)) + geom_histogram(fill="grey88", col="black",breaks=seq(0, 1, by = 0.1)) + xlim(0,1) + ylim(0,23) + theme(axis.text=element_blank(), axis.title=element_blank(), axis.line=element_blank(), axis.ticks=element_blank())
points_vaujany<-ggplot(table_vaujany, aes(HINDEX, Heterozygosity)) + geom_point(size=6, position="jitter") + xlim(0,1) + ylim(0,0.5) + background_grid(major = 'y', minor = "none") + panel_border() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())
hist_vaujany<-ggplot(table_vaujany, aes(HINDEX)) + geom_histogram(fill="grey88", col="black",breaks=seq(0, 1, by = 0.1)) + xlim(0,1) + ylim(0,23) + theme(axis.text=element_blank(), axis.title=element_blank(), axis.line=element_blank(), axis.ticks=element_blank())
points_locarno<-ggplot(table_locarno, aes(HINDEX, Heterozygosity)) + geom_point(size=6, position="jitter") + xlim(0,1) + ylim(0,0.5) + background_grid(major = 'y', minor = "none") + panel_border() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())
hist_locarno<-ggplot(table_locarno, aes(HINDEX)) + geom_histogram(fill="grey88", col="black",breaks=seq(0, 1, by = 0.1)) + xlim(0,1) + ylim(0,23) + theme(axis.text=element_blank(), axis.title=element_blank(), axis.line=element_blank(), axis.ticks=element_blank())
points_chaillol<-ggplot(table_chaillol, aes(HINDEX, Heterozygosity)) + geom_point(size=6, position="jitter") + xlim(0,1) + ylim(0,0.5) + background_grid(major = 'y', minor = "none") + panel_border() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())
hist_chaillol<-ggplot(table_chaillol, aes(HINDEX)) + geom_histogram(fill="grey88", col="black",breaks=seq(0, 1, by = 0.1)) + xlim(0,1) + ylim(0,23) + theme(axis.text=element_blank(), axis.title=element_blank(), axis.line=element_blank(), axis.ticks=element_blank())
points_vars<-ggplot(table, aes(HINDEX, Heterozygosity)) + geom_point(size=6, position="jitter") + xlim(0,1) + ylim(0,0.5) + background_grid(major = 'y', minor = "none") + panel_border() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())
hist_vars<-ggplot(table, aes(HINDEX)) + geom_histogram(fill="grey88", col="black",breaks=seq(0, 1, by = 0.1)) + xlim(0,1) + ylim(0,23) + theme(axis.text=element_blank(), axis.title=element_blank(), axis.line=element_blank(), axis.ticks=element_blank())
points_bel<-ggplot(table_bel, aes(HINDEX, Heterozygosity)) + geom_point(size=6, position="jitter") + xlim(0,1) + ylim(0,0.5) + background_grid(major = 'y', minor = "none") + panel_border() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())
hist_bel<-ggplot(table_bel, aes(HINDEX)) + geom_histogram(fill="grey88", col="black",breaks=seq(0, 1, by = 0.1)) + xlim(0,1) + ylim(0,23) + theme(axis.text=element_blank(), axis.title=element_blank(), axis.line=element_blank(), axis.ticks=element_blank())

points_vaujany2<-ggplot(table_vaujany, aes(HINDEX, Heterozygosity)) + geom_point(size=6, position="jitter") + xlim(0,1) + ylim(0,0.5) + background_grid(major = 'y', minor = "none") + panel_border()
hist_vaujany2<-ggplot(table_vaujany, aes(HINDEX)) + geom_histogram(fill="grey88", col="black",breaks=seq(0, 1, by = 0.1)) + xlim(0,1) + ylim(0,23) + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank())

pdf("HINDEX_Heterozigosity.pdf", width=20, height=5)
grid.arrange(hist_vaujany2, hist_vaujany, hist_aoste, hist_chaillol, hist_locarno, hist_vars, hist_bel, points_vaujany2, points_vaujany, points_aoste, points_chaillol, points_locarno, points_vars, points_bel, ncol=7,nrow=2, heights=c(2, 2), widths=rep(1.4,7))
dev.off()



toto<-ggplot(table, aes(HINDEX, Heterozygosity)) + geom_point(size=5, position="jitter") + xlim(0,1) + ylim(0,0.5) + theme_bw(30)
ggExtra::ggMarginal(toto, type = "histogram", margin = 'x', size=1.5, xparams = list(binwidth = 0.1, fill = "grey88"))




#### Comparaison des deux clines ####

compare.clines(cline_out_ecrins,cline_out_vars,sig.test=TRUE,n.reps=1000)

multinom(cline_out_ecrins$Count.matrix[1,]~cline_out_ecrins$hybrid.index)

##########################################################
#### Représentation graphique HINDEX et heterozygotie ####

## Distribution des H index sur les differentes pops etudiees

par(mfrow=c(3,2))
hist(HINDEX_vars$h, breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
hist(HINDEX_ecrins$h, breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
hist(HINDEX_vaujany$h, breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
hist(HINDEX_aoste$h, breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
hist(HINDEX_chaillol$h, breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")

pdf(paste("C:/Users/capblancq/Documents/Projet Papillons/These Hybridation Coenonympha/ddRADseq/Analyses ddRADseq/HINDEX/Hist_Hindex_Vaujany", ".pdf",sep=""))
hist(HINDEX_vaujany$h, breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
dev.off()


## Histogramme de l'hétérozygotie dans chaque pop
Loc_vars<-Loc[which(Loc=="SLC"|Loc=="MEL"|Loc=="SPA"|Loc=="VAR")] ## Vecteur populations correspondant aux indice d'heterozygotie
barplot(c(mean(het_vars[which(Loc_vars=="VAR")]),mean(het_vars[which(Loc_vars=="SLC")]),mean(het_vars[which(Loc_vars=="MEL")]),mean(het_vars[which(Loc_vars=="SPA")]))) ## Representation graphique

Loc_ecrins<-Loc[which(Loc=="AIL"|Loc=="NAV"|Loc=="GIO"|Loc=="CHL_dar"|Loc=="ORC_dar"|Loc=="BY"|Loc=="VALS"|Loc=="SEL"|Loc=="DOR")] ## Vecteur populations correspondant aux indice d'heterozygotie
barplot(c(mean(het_ecrins[which(Loc_ecrins=="AIL")]),mean(het_ecrins[which(Loc_ecrins=="SEL")]),mean(het_ecrins[which(Loc_ecrins=="BY")]),mean(het_ecrins[which(Loc_ecrins=="VALS")]),mean(het_ecrins[which(Loc_ecrins=="NAV")]),mean(het_ecrins[which(Loc_ecrins=="GIO")]),mean(het_ecrins[which(Loc_ecrins=="DOR")]),mean(het_ecrins[which(Loc_ecrins=="ORC_dar")]),mean(het_ecrins[which(Loc_ecrins=="CHL_dar")]))) ## Representation graphique

## Histogramme des H INDEX par pop sur le massif des Ecrins

par(mfrow=c(3,3))
hist(HINDEX_ecrins$h[which(Loc_ecrins=="AIL")], breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey", axes=F)
hist(HINDEX_ecrins$h[which(Loc_ecrins=="NAV")], breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
hist(HINDEX_ecrins$h[which(Loc_ecrins=="GIO")], breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
hist(HINDEX_ecrins$h[which(Loc_ecrins=="CHL_dar")], breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
hist(HINDEX_ecrins$h[which(Loc_ecrins=="ORC_dar")], breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
hist(HINDEX_ecrins$h[which(Loc_ecrins=="BY")], breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
hist(HINDEX_ecrins$h[which(Loc_ecrins=="VALS")], breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
hist(HINDEX_ecrins$h[which(Loc_ecrins=="DOR")], breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
hist(HINDEX_ecrins$h[which(Loc_ecrins=="SEL")], breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")

## Histogramme des H INDEX par pop sur le col de Vars

par(mfrow=c(1,4))
hist(HINDEX_vars$h[which(Loc_vars=="VAR")], breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
hist(HINDEX_vars$h[which(Loc_vars=="SLC")], breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
hist(HINDEX_vars$h[which(Loc_vars=="MEL")], breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")
hist(HINDEX_vars$h[which(Loc_vars=="SPA")], breaks=seq(0,1,0.1), freq=T, xlim=c(0,1), col="darkgrey")


triangle_plot_bis<-function (hi.index = NULL, int.het = NULL, pdf = TRUE, out.file = "tri_plot.pdf") 
{
  if (is.null(hi.index) == TRUE | is.null(int.het) == TRUE) 
    stop("error, input data are not provided")
  if (is.data.frame(hi.index) == TRUE) 
    hi.index <- hi.index[, 2]
  if (pdf == TRUE) 
    pdf(file = paste(out.file))
  plot(hi.index, int.het, xlab = "Hybrid index", ylab = "Interspecific heterozygosity", 
       xlim = c(0, 1), ylim = c(0, 1), cex=4, col="grey30", pch=16)
  lines(c(0, 0.5), c(0, 1))
  lines(c(0.5, 1), c(1, 0))
  if (pdf == TRUE) 
    dev.off()
}

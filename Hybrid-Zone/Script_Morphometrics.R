##############################################
#### Geometric Morphometrics Hybrid zones ####
##############################################

setwd("./")

## Packages
library(ade4)
library(MASS)
library(ape)
library(ggplot2)

## Functions of Rmorpho package
source("./Fonctions du package Rmorph/evolCVP_multi.R")
source("./Fonctions du package Rmorph/FonctionsR-Alloselection.R")
source("./Fonctions du package Rmorph/FonctionsR-Alloselection.R2.R")
source("./Fonctions du package Rmorph/Rmorph_Rapplic_release_v04-12.R")
source("./Fonctions du package Rmorph/Rmorph_devel_Rapplic_2010.R")

##############
#### DATA ####

## Chargement des fichiers de coordonn?es pour chaque aile et des informations sur les individus
coordAP<-read.table("Coord_PostD_all.txt")
infoAP<-read.table("Ind_PostD_all.txt")
datAP<-col2mat(coordAP,npts=22)
tableAP<-cbind(infoAP,datAP)
coordVA<-read.table("Coord_AntD_all.txt")
infoVA<-read.table("Ind_AntD_all.txt")
datVA<-col2mat(coordVA,npts=18)
tableVA<-cbind(infoVA,datVA)

## Merging hind and fore-wing tables
table_all<-merge(tableAP,tableVA,by="V1")
table_all<-table_all[,-which(colnames(table_all)=="V2.y")]
colnames(table_all)<-c("Id", "Pop",as.character(1:80))

## Populations related to the hybrid zone
info<-read.table("./infos_ind_VARS", header=T, sep = ",")
table_Vars <- table_all[which(table_all$Pop%in%info$Pop),]
table_Vars <- rbind(table_Vars[which(table_Vars$Pop=="AOS"|table_Vars$Pop=="LAU"|table_Vars$Pop=="SES"|table_Vars$Pop=="AIL"|table_Vars$Pop=="VJ_7"|table_Vars$Pop=="VJ_8"|table_Vars$Pop=="VJ_9"),], table_Vars[which(table_Vars$Pop=="LAR"|table_Vars$Pop=="BOR"|table_Vars$Pop=="FOA"|table_Vars$Pop=="LOM"|table_Vars$Pop=="SEY"),], table_Vars[which(table_Vars$Pop=="TOU"|table_Vars$Pop=="SLC"|table_Vars$Pop=="MEL_INF"|table_Vars$Pop=="MEL_SUP"|table_Vars$Pop=="SPA_1"|table_Vars$Pop=="SPA_2"|table_Vars$Pop=="SPA_3"|table_Vars$Pop=="VAR_1"|table_Vars$Pop=="VAR_2"|table_Vars$Pop=="VAR_3"|table_Vars$Pop=="VAR_N"),])
table_Vars <- table_Vars[-which(table_Vars$Id%in%c("AIL_10", "SEY_18", "AOS_1", "AOS_2", "AOS_3", "AOS_4", "AOS_5","AOS_12","AOS_13","AOS_15","AOS_18","AOS_20","AOS_21")),] ## remove C. arcania individuals

#####################################################
#### Splitting eyespots, venation and white band ####

supAP_ven<-gpa(table_Vars[,c(3:16,35:38)], dimx=2)
scoresAP_ven<-supAP_ven$scores
supAP_wb<-gpa(table_Vars[,17:34], dimx=2)
scoresAP_wb<-supAP_wb$scores
supAP_eye<-gpa(table_Vars[,39:46], dimx=2)
scoresAP_eye<-supAP_eye$scores
supAA_ven<-gpa(table_Vars[,47:82], dimx=2)
scoresAA_ven<-supAA_ven$scores

## Differences of shape
visu(supAP_ven, scoresAP_ven[,1], GP = "BW")
visu(supAA_ven, scoresAA_ven[,1], GP = "BW")
visu(supAP_wb, scoresAP_wb[,1], GP = "BW")
visu(supAP_eye, scoresAP_eye[,1], GP = "BW")

## Sorting by populations
scoresAP_ven_parents<-scoresAP_ven[which(table_Vars$Pop=="AOS"|table_Vars$Pop=="LAU"|table_Vars$Pop=="SES"|table_Vars$Pop=="AIL"|table_Vars$Pop=="VJ_7"|table_Vars$Pop=="VJ_8"|table_Vars$Pop=="VJ_9"|table_Vars$Pop=="LAR"|table_Vars$Pop=="BOR"|table_Vars$Pop=="FOA"|table_Vars$Pop=="LOM"|table_Vars$Pop=="SEY"),]
scoresAA_ven_parents<-scoresAA_ven[which(table_Vars$Pop=="AOS"|table_Vars$Pop=="LAU"|table_Vars$Pop=="SES"|table_Vars$Pop=="AIL"|table_Vars$Pop=="VJ_7"|table_Vars$Pop=="VJ_8"|table_Vars$Pop=="VJ_9"|table_Vars$Pop=="LAR"|table_Vars$Pop=="BOR"|table_Vars$Pop=="FOA"|table_Vars$Pop=="LOM"|table_Vars$Pop=="SEY"),]
scoresAP_eye_parents<-scoresAP_eye[which(table_Vars$Pop=="AOS"|table_Vars$Pop=="LAU"|table_Vars$Pop=="SES"|table_Vars$Pop=="AIL"|table_Vars$Pop=="VJ_7"|table_Vars$Pop=="VJ_8"|table_Vars$Pop=="VJ_9"|table_Vars$Pop=="LAR"|table_Vars$Pop=="BOR"|table_Vars$Pop=="FOA"|table_Vars$Pop=="LOM"|table_Vars$Pop=="SEY"),]
scoresAP_wb_parents<-scoresAP_wb[which(table_Vars$Pop=="AOS"|table_Vars$Pop=="LAU"|table_Vars$Pop=="SES"|table_Vars$Pop=="AIL"|table_Vars$Pop=="VJ_7"|table_Vars$Pop=="VJ_8"|table_Vars$Pop=="VJ_9"|table_Vars$Pop=="LAR"|table_Vars$Pop=="BOR"|table_Vars$Pop=="FOA"|table_Vars$Pop=="LOM"|table_Vars$Pop=="SEY"),]

scoresAP_ven_hybrids<-scoresAP_ven[which(table_Vars$Pop=="TOU"|table_Vars$Pop=="SLC"|table_Vars$Pop=="MEL_INF"|table_Vars$Pop=="MEL_SUP"|table_Vars$Pop=="SPA_1"|table_Vars$Pop=="SPA_2"|table_Vars$Pop=="SPA_3"|table_Vars$Pop=="VAR_1"|table_Vars$Pop=="VAR_2"|table_Vars$Pop=="VAR_3"|table_Vars$Pop=="VAR_N"),]
scoresAA_ven_hybrids<-scoresAA_ven[which(table_Vars$Pop=="TOU"|table_Vars$Pop=="SLC"|table_Vars$Pop=="MEL_INF"|table_Vars$Pop=="MEL_SUP"|table_Vars$Pop=="SPA_1"|table_Vars$Pop=="SPA_2"|table_Vars$Pop=="SPA_3"|table_Vars$Pop=="VAR_1"|table_Vars$Pop=="VAR_2"|table_Vars$Pop=="VAR_3"|table_Vars$Pop=="VAR_N"),]
scoresAP_eye_hybrids<-scoresAP_eye[which(table_Vars$Pop=="TOU"|table_Vars$Pop=="SLC"|table_Vars$Pop=="MEL_INF"|table_Vars$Pop=="MEL_SUP"|table_Vars$Pop=="SPA_1"|table_Vars$Pop=="SPA_2"|table_Vars$Pop=="SPA_3"|table_Vars$Pop=="VAR_1"|table_Vars$Pop=="VAR_2"|table_Vars$Pop=="VAR_3"|table_Vars$Pop=="VAR_N"),]
scoresAP_wb_hybrids<-scoresAP_wb[which(table_Vars$Pop=="TOU"|table_Vars$Pop=="SLC"|table_Vars$Pop=="MEL_INF"|table_Vars$Pop=="MEL_SUP"|table_Vars$Pop=="SPA_1"|table_Vars$Pop=="SPA_2"|table_Vars$Pop=="SPA_3"|table_Vars$Pop=="VAR_1"|table_Vars$Pop=="VAR_2"|table_Vars$Pop=="VAR_3"|table_Vars$Pop=="VAR_N"),]

par(mfrow=c(2,3))
plotg(scoresAP_ven_parents, c(1,2), group=Group_parents, GP="col", colpoints=c("#3288BD","#FEE08B"),labX="PC1 (32%)", labY="PC2 (20%)", gtitle="Hindwing venation") 
plotg(scoresAA_ven_parents, c(1,2), group=Group_parents, GP="col", colpoints=c("#3288BD","#FEE08B"),labX="PC1 (35%)", labY="PC2 (14%)", gtitle="Forewing venation")
plotg(scoresAP_eye_parents, c(1,2), group=Group_parents, GP="col", colpoints=c("#3288BD","#FEE08B"),labX="PC1 (66%)", labY="PC2 (18%)", gtitle="Eyespots alignment") 
plotg(scoresAP_wb_parents, c(1,2), group=Group_parents, GP="col", colpoints=c("#3288BD","#FEE08B"),labX="PC1 (47%)", labY="PC2 (26%)", gtitle="White band shape") 
boxplot((CSizeVars[1:length(Group_parents),1]+CSizeVars[1:length(Group_parents),2])~as.character(Group_parents), col=c("#FEE08B", "#3288BD"), main = "Relative wings size")
dev.off()


summary(manova(scoresAP_ven_parents[,1:2]~Group_parents))
summary(manova(scoresAA_ven_parents[,1:2]~Group_parents))
summary(manova(scoresAP_eye_parents[,1:2]~Group_parents))
summary(manova(scoresAP_wb_parents[,1:2]~Group_parents))

t.test((CSizeVars[1:length(Group_parents),1]+CSizeVars[1:length(Group_parents),2])~as.character(Group_parents))$p.value

#boxplot(as.numeric(scoresAP_eye_hybrids[,1])~as.factor(as.character(Group_hybrids)))

##################################
#### TABLE for cline analysis ####

TAB <- cbind(table_Vars[,1:2], scoresAP_ven[,1], scoresAA_ven[,1], scoresAP_eye[,1], scoresAP_wb[,1], CSizeVars, WingSize = (CSizeVars[,1]+CSizeVars[,2]))
TAB_pop_mean <- aggregate(TAB[,-c(1:2)], by=list(TAB$Pop), FUN="mean", na.rm=TRUE)
TAB_pop_var <- aggregate(TAB[,-c(1:2)], by=list(TAB$Pop), FUN="var", na.rm=TRUE)
TAB_pop_count <- aggregate(TAB$Pop, by=list(TAB$Pop), FUN="length")

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

TAB_scaled <- cbind(table_Vars[,1:2], Eyespots = range01(scoresAP_eye[,1]), Whiteband = range01(scoresAP_wb[,1]), WingSize = range01(CSizeVars[,1]+CSizeVars[,2]))
TAB_pop_scaled_mean <- aggregate(TAB_scaled[,-c(1:2)], by=list(TAB_scaled$Pop), FUN="mean", na.rm=TRUE)
TAB_pop_scaled_var <- aggregate(TAB_scaled[,-c(1:2)], by=list(TAB_scaled$Pop), FUN="var", na.rm=TRUE)
TAB_pop_scaled_count <- aggregate(TAB_scaled$Pop, by=list(TAB_scaled$Pop), FUN="length")

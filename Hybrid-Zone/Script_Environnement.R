##############################################################################################
##############################################################################################

# Script genetic and environment association (Coenonypmha) 
# Capblancq Thibaut
# Date : 30-10-2018

##############################################################################################
##############################################################################################

library(dismo)
library(ade4)
library(raster)
library(SDMTools)
library(sm)
library(ggplot2)
library(reshape2)

setwd("./")

##############
#### DATA ####

## Infos
info<-read.table("./infos_ind_VARS.txt", header=T, sep = ",")
  
## rasters environment
fich_var_1<-list.files(path='./Env_variables/',pattern='FA.tif$',full.names=T)
fich_var_2<-list.files(path='./Env_variables/',pattern='.img$',full.names=T)
ras_env <- stack(fich_var_1)
ras_env <- stack(ras_env, projectRaster(stack(fich_var_2), ras_env))

## Env data extraction
var_env <- extract(ras_env,info[,4:5])
TAB <- cbind(info[,1:3], var_env)
TAB_pop <- aggregate(TAB, list(info$Pop), FUN = head, 1)

## HINDEX
TAB$HINDEX <- rep(NA, length(TAB$Ind))
TAB$HINDEX[which(TAB$Zone=="GARDETTA")] <- 1
TAB$HINDEX[which(TAB$Zone=="MACROMMA")] <- 0
TAB$HINDEX[which(TAB$Zone=="VARS" & TAB$Ind %in% row.names(HINDEX_VARS))] <- HINDEX_VARS$h[match(TAB$Ind[which(TAB$Zone=="VARS" & TAB$Ind %in% row.names(HINDEX_VARS))], row.names(HINDEX_VARS))]
TAB <- TAB[-is.na(TAB$HINDEX),]

## Parental species divergence
TAB_temp <- TAB[which(TAB$Zone=="GARDETTA" | TAB$Zone=="MACROMMA"),]
pdf("./Diff_vareco.pdf", width=10, height=6)
par(mfrow=c(2,4))
p.val <- NULL
for(i in 4:11)
{
  boxplot(TAB_temp[,i]~as.character(TAB_temp$Zone), col=c("#3288BD", "#FEE08B"), main=colnames(TAB_temp)[i])
  p.val <- c(p.val, t.test(TAB_temp[,i]~as.character(TAB_temp$Zone))$p.value)
}
dev.off()


##############
#### Vars ####

TAB_VARS <- TAB[which(TAB$Zone=="VARS"),c(1,3,4,5,7,8,9,10,11)]
info_VARS <- info[which(info$Zone=="VARS"),]
TAB_VARS <- TAB_VARS[order(info_VARS$Y_Lambert93),]

## biovolume
biovolume_VARS <- read.table("./Env_variables/Variables_Biovolume_Vars.txt", header=T)
biovolume_VARS <- aggregate(biovolume_VARS, list(biovolume_VARS$alt), mean)[,c("NB_TrSh", "H_max", "biovolume", "herb", "grass")]
biovolume_VARS[11,] <- c(0, NA, NA, NA, NA)

## Distance along the hybrid zone coming from Script_Genomic_Clines.R
dist_Vars<-pointDistance(c(997027.1, 6383629.52),info_VARS[,4:5], lonlat=FALSE)
dist_Vars_pop <- aggregate(dist_Vars, by = list(info_VARS$Pop), mean)
dist_Vars_pop <- dist_Vars_pop[order(dist_Vars_pop$x),]

TAB_VARS$Nb_TreeShrubs <- rep(biovolume_VARS$NB_TrSh, times=as.vector(aggregate(TAB_VARS$Pop, by = list(TAB_VARS$Pop), FUN=length)[,2]))

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

TAB_VARS <- data.frame(TAB_VARS[,1:2], apply(TAB_VARS[,-c(1:2)], 2, range01))

#########################
#### Cline with hzar ####

registerDoParallel(cores = 8)
c.fit.env = foreach (i = colnames(TAB_VARS[,-c(1:2)])) %dopar%
{

c.env <- list()
## Space to hold the observed data
c.env$obs <- list()
## Space to hold the models to fit
c.env$models <- list()
## Space to hold the compiled fit requests
c.env$fitRs <- list()
## Space to hold the output data chains
c.env$runs <- list()
## Space to hold the analysed data
c.env$analysis <- list()

## Formating data
c.env$obs <- hzar.doNormalData1DRaw(hzar.mapSiteDist(dist_Vars_pop$Group.1,dist_Vars_pop$x),TAB_VARS$Pop,TAB_VARS[,i])

## Looking at the observed data
# hzar.plot.obsData(c.env$obs)

## Creating models 
c.env.loadmodel <- function(scaling,tails,id=paste(scaling,tails,sep="."))
{
  c.env$models[[id]] <<- hzar.makeCline1DNormal(c.env$obs, tails)
  if (all(regexpr("fixed",scaling,ignore.case=TRUE) == 1 ))
  {hzar.meta.fix(c.env$models[[id]])$muL <<- TRUE
  hzar.meta.fix(c.env$models[[id]])$muR <<- TRUE
  hzar.meta.fix(c.env$models[[id]])$varL <<- TRUE
  hzar.meta.fix(c.env$models[[id]])$varR <<- TRUE
  }
  hzar.meta.init(c.env$models[[id]])$muL <<- c.env$obs$frame["TOU","mu"]
  #hzar.meta.init(c.env$models[[id]])$varL <<- c.env$obs$frame["TOU","var"]
  hzar.meta.init(c.env$models[[id]])$muR <<- c.env$obs$frame["VAR_N","mu"]
  #hzar.meta.init(c.env$models[[id]])$varR <<- c.env$obs$frame["VAR_N","var"]
}

c.env.loadmodel("fixed","none","model.none")
c.env.loadmodel("fixed","both","model.both")
c.env.loadmodel("fixed","right","model.right")
c.env.loadmodel("fixed","left","model.left")
c.env.loadmodel("fixed","mirror","model.mirror")
c.env.loadmodel("free","none","modelII")
c.env.loadmodel("free","both","modelIII")

## Modifying model parameters
# Bounds
c.env$models <- sapply(c.env$models,hzar.model.addBoxReq, 0, 9000, simplify=FALSE)
# reduce the tune setting of modelII from 1.5 to 1.4
hzar.meta.tune(c.env$models$modelII)<-1.4
# reduce the tune setting of modelIII from 1.5 to 1.2
hzar.meta.tune(c.env$models$modelIII)<-1.2

## Compile models to prepare for fitting
c.env$fitRs$init <- sapply(c.env$models, hzar.first.fitRequest.gC, obsData=c.env$obs, verbose=FALSE, simplify=FALSE)

## Model chain length and seeds
chainLength=1e5;                       
mainSeed=list(A=c(978,544,99,596,528,124), B=c(544,99,596,528,124,978), C=c(99,596,528,124,978,544), D=c(978,99,596,528,124,544), E=c(124,596,978,99,528,544), F=c(99,124,596,978,528,544), G=c(528,99,124,596,978,544))
lapply(c.env$fitRs$init, function(x) x$mcmcParam$chainLength <- chainLength)
lapply(c.env$fitRs$init, function(x) x$mcmcParam$burnin <- chainLength %/% 10)
lapply(1:length(c.env$fitRs$init), function(x) c.env$fitRs$init[[x]]$mcmcParam$seed[[1]] <- mainSeed[[x]])

## Run model for an initial chain
c.env$runs$init <- list()
c.env$runs$init$model.none <- hzar.doFit(c.env$fitRs$init$model.none)
c.env$runs$init$model.both <- hzar.doFit(c.env$fitRs$init$model.both)
c.env$runs$init$model.right <- hzar.doFit(c.env$fitRs$init$model.right)
c.env$runs$init$model.left <- hzar.doFit(c.env$fitRs$init$model.left)
c.env$runs$init$model.mirror <- hzar.doFit(c.env$fitRs$init$model.mirror)
c.env$runs$init$modelII <- hzar.doFit(c.env$fitRs$init$modelII)
c.env$runs$init$modelIII <- hzar.doFit(c.env$fitRs$init$modelIII)

## Compile a new set of fi requests using the initial chains
c.env$fitRs$chains <- lapply(c.env$runs$init, hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
c.env$fitRs$chains <- hzar.multiFitRequest(c.env$fitRs$chains, each=3, baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit
# center for all models
lapply(1:length(c.env$fitRs$chains), function(x) c.env$fitRs$chains[[x]]$modelParam$init["center"]=runif(21,0,9000)[x])
# width for all models
lapply(1:length(c.env$fitRs$chains), function(x) c.env$fitRs$chains[[x]]$modelParam$init["width"]=runif(21,0,10000)[x])
# varH for all models
lapply(1:length(c.env$fitRs$chains), function(x) c.env$fitRs$chains[[x]]$modelParam$init["varH"]=10^runif(21,-1,1)[x])

## Go ahead and run a chain of 3 runs for every fit request
c.env$runs$chains <-  hzar.doChain.multi(c.env$fitRs$chains, doPar=TRUE, inOrder=FALSE, count=3)

## Start aggregation of data for analysis
c.env$analysis$initDGs <- list()

## Create a model data group for each model
c.env$analysis$initDGs$model.none <- hzar.dataGroup.add(c.env$runs$init$model.none)
c.env$analysis$initDGs$model.both <- hzar.dataGroup.add(c.env$runs$init$model.both)
c.env$analysis$initDGs$model.right <- hzar.dataGroup.add(c.env$runs$init$model.right)
c.env$analysis$initDGs$model.left<- hzar.dataGroup.add(c.env$runs$init$model.left)
c.env$analysis$initDGs$model.mirror <- hzar.dataGroup.add(c.env$runs$init$model.mirror)
c.env$analysis$initDGs$modelII <- hzar.dataGroup.add(c.env$runs$init$modelII)
c.env$analysis$initDGs$modelIII <- hzar.dataGroup.add(c.env$runs$init$modelIII)

## Create a hzar.obsDataGroup object from the models
c.env$analysis$oDG <- hzar.make.obsDataGroup(c.env$analysis$initDGs)
c.env$analysis$oDG <- hzar.copyModelLabels(c.env$analysis$initDGs, c.env$analysis$oDG)

## Convert all runs to hzar.dataGroup objects
c.env$analysis$oDG <- hzar.make.obsDataGroup(lapply(c.env$runs$chains, hzar.dataGroup.add), c.env$analysis$oDG)

## Do model selection based on the AICs scores
c.env$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(c.env$analysis$oDG)

## Print out the model with the minimum AICc score
c.env$analysis$model.name <- rownames(c.env$analysis$AICcTable)[[which.min(c.env$analysis$AICcTable$AICc)]]

## Extract the hzar.dataGroup object for the selected model
c.env$analysis$model.selected <- c.env$analysis$oDG$data.groups[[c.env$analysis$model.name]]

return(c.env)

}

########################
#### Model selected ####

models <- as.data.frame(lapply(c.fit.env, function(x) c(unlist(which.min(hzar.AICc.hzar.obsDataGroup(x$analysis$oDG)[,1])), unlist(hzar.AICc.hzar.obsDataGroup(x$analysis$oDG)[which.min(hzar.AICc.hzar.obsDataGroup(x$analysis$oDG)[,1]),]))))
colnames(models) <- colnames(TAB_VARS[,-c(1:2)])


##########################
#### cline parameters ####
c<-NULL
w<-NULL
width<-NULL
center<-NULL
slope<-NULL
par<-NULL
params_VARS<-NULL
for (i in 1:length(c.fit.env))
{
  c<-density(c.fit.env[[i]]$analysis$model.selected$data.param$center)
  center<-c$x[which.max(c$y)]
  w<-density(c.fit.env[[i]]$analysis$model.selected$data.param$width)
  width<-w$x[which.max(w$y)]
  slope<-(abs(c.fit.env[[i]]$analysis$model.selected$ML.cline$param.all$muR-c.fit.env[[i]]$analysis$model.selected$ML.cline$param.all$muL)/width)
  par<-c(center, width, slope)
  params_VARS<-rbind(params_VARS, par)
}
row.names(params_VARS)<-colnames(TAB_VARS[,-c(1:2)])
colnames(params_VARS)<-c("center", "width", "Slope")
params_VARS<-as.data.frame(params_VARS)
write.table(params_VARS, "params_clines_env.csv", sep=",")

### Visualization
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

## Changing cline bounds to compare them
cline.env.scaled <- list()
cline.env.scaled[[1]] <- hzar.gen.cline(list(center = c.fit.env[[1]]$analysis$model.selected$ML.cline$param.all$center, width = c.fit.env[[1]]$analysis$model.selected$ML.cline$param.all$width, muL=0,muR=1, varH = c.fit.env[[1]]$analysis$model.selected$ML.cline$param.all$varH, varL = c.fit.env[[1]]$analysis$model.selected$ML.cline$param.all$varL, varR = c.fit.env[[1]]$analysis$model.selected$ML.cline$param.all$varR), c.fit.env[[1]]$analysis$model.selected)
cline.env.scaled[[2]] <- hzar.gen.cline(list(center = c.fit.env[[2]]$analysis$model.selected$ML.cline$param.all$center, width = c.fit.env[[2]]$analysis$model.selected$ML.cline$param.all$width, deltaL = c.fit.env[[2]]$analysis$model.selected$ML.cline$param.all$deltaL, tauL = c.fit.env[[2]]$analysis$model.selected$ML.cline$param.all$tauL, deltaR = c.fit.env[[2]]$analysis$model.selected$ML.cline$param.all$deltaR, tauR = c.fit.env[[2]]$analysis$model.selected$ML.cline$param.all$tauR, muL=1,muR=0, varH = c.fit.env[[2]]$analysis$model.selected$ML.cline$param.all$varH, varL = c.fit.env[[2]]$analysis$model.selected$ML.cline$param.all$varL, varR = c.fit.env[[2]]$analysis$model.selected$ML.cline$param.all$varR), c.fit.env[[2]]$analysis$model.selected)
cline.env.scaled[[3]] <- hzar.gen.cline(list(center = c.fit.env[[3]]$analysis$model.selected$ML.cline$param.all$center, width = c.fit.env[[3]]$analysis$model.selected$ML.cline$param.all$width, deltaL = c.fit.env[[3]]$analysis$model.selected$ML.cline$param.all$deltaL, tauL = c.fit.env[[3]]$analysis$model.selected$ML.cline$param.all$tauL, deltaR = c.fit.env[[3]]$analysis$model.selected$ML.cline$param.all$deltaR, tauR = c.fit.env[[3]]$analysis$model.selected$ML.cline$param.all$tauR, muL=1,muR=0, varH = c.fit.env[[3]]$analysis$model.selected$ML.cline$param.all$varH, varL = c.fit.env[[3]]$analysis$model.selected$ML.cline$param.all$varL, varR = c.fit.env[[3]]$analysis$model.selected$ML.cline$param.all$varR), c.fit.env[[3]]$analysis$model.selected)
cline.env.scaled[[4]] <- hzar.gen.cline(list(center = c.fit.env[[4]]$analysis$model.selected$ML.cline$param.all$center, width = c.fit.env[[4]]$analysis$model.selected$ML.cline$param.all$width, muL=0,muR=1, varH = c.fit.env[[4]]$analysis$model.selected$ML.cline$param.all$varH, varL = c.fit.env[[4]]$analysis$model.selected$ML.cline$param.all$varL, varR = c.fit.env[[4]]$analysis$model.selected$ML.cline$param.all$varR), c.fit.env[[4]]$analysis$model.selected)
cline.env.scaled[[5]] <- hzar.gen.cline(list(center = c.fit.env[[5]]$analysis$model.selected$ML.cline$param.all$center, width = c.fit.env[[5]]$analysis$model.selected$ML.cline$param.all$width, deltaL = c.fit.env[[5]]$analysis$model.selected$ML.cline$param.all$deltaL, tauL = c.fit.env[[5]]$analysis$model.selected$ML.cline$param.all$tauL, deltaR = c.fit.env[[5]]$analysis$model.selected$ML.cline$param.all$deltaR, tauR = c.fit.env[[5]]$analysis$model.selected$ML.cline$param.all$tauR, muL=1,muR=0, varH = c.fit.env[[5]]$analysis$model.selected$ML.cline$param.all$varH, varL = c.fit.env[[5]]$analysis$model.selected$ML.cline$param.all$varL, varR = c.fit.env[[5]]$analysis$model.selected$ML.cline$param.all$varR), c.fit.env[[5]]$analysis$model.selected)
cline.env.scaled[[6]] <- c.fit.env[[6]]$analysis$model.selected #hzar.gen.cline(list(center = c.fit.env[[6]]$analysis$model.selected$ML.cline$param.all$center, width = c.fit.env[[6]]$analysis$model.selected$ML.cline$param.all$width, deltaL = c.fit.env[[6]]$analysis$model.selected$ML.cline$param.all$deltaL, tauL = c.fit.env[[6]]$analysis$model.selected$ML.cline$param.all$tauL, deltaR = c.fit.env[[6]]$analysis$model.selected$ML.cline$param.all$deltaR, tauR = c.fit.env[[6]]$analysis$model.selected$ML.cline$param.all$tauR, muL=1,muR=0, varH = c.fit.env[[6]]$analysis$model.selected$ML.cline$param.all$varH, varL = c.fit.env[[6]]$analysis$model.selected$ML.cline$param.all$varL, varR = c.fit.env[[6]]$analysis$model.selected$ML.cline$param.all$varR), c.fit.env[[6]]$analysis$model.selected)
cline.env.scaled[[7]] <- hzar.gen.cline(list(center = c.fit.env[[7]]$analysis$model.selected$ML.cline$param.all$center, width = c.fit.env[[7]]$analysis$model.selected$ML.cline$param.all$width, deltaL = c.fit.env[[7]]$analysis$model.selected$ML.cline$param.all$deltaL, tauL = c.fit.env[[7]]$analysis$model.selected$ML.cline$param.all$tauL, deltaR = c.fit.env[[7]]$analysis$model.selected$ML.cline$param.all$deltaR, tauR = c.fit.env[[7]]$analysis$model.selected$ML.cline$param.all$tauR, muL=1,muR=0, varH = c.fit.env[[7]]$analysis$model.selected$ML.cline$param.all$varH, varL = c.fit.env[[7]]$analysis$model.selected$ML.cline$param.all$varL, varR = c.fit.env[[7]]$analysis$model.selected$ML.cline$param.all$varR), c.fit.env[[7]]$analysis$model.selected)
cline.env.scaled[[8]] <- hzar.gen.cline(list(center = c.fit.env[[8]]$analysis$model.selected$ML.cline$param.all$center, width = c.fit.env[[8]]$analysis$model.selected$ML.cline$param.all$width, muL=1,muR=0, varH = c.fit.env[[8]]$analysis$model.selected$ML.cline$param.all$varH, varL = c.fit.env[[8]]$analysis$model.selected$ML.cline$param.all$varL, varR = c.fit.env[[8]]$analysis$model.selected$ML.cline$param.all$varR), c.fit.env[[8]]$analysis$model.selected)

cols <- c(brewer.pal(8, "Dark2")[-c(1,4)], "#3288BD","#66C2A5")

pdf("env_clines_Vars_bounded.pdf", width = 5, height = 5)
plot(1,1, xlim=c(0,8000), ylim=c(0,1), "n", ylab = "Percent characteristic value", xlab = "Distance (m)")
for(i in c(1:8)){hzar.plot.cline2(cline.env.scaled[[i]], col=cols[i], add=TRUE, lwd=4, lty=1)}
legend(6000,1, legend = colnames(TAB_VARS[,-c(1:2)]), col = cols, lty=1, lwd=4, cex = 0.8)
dev.off()

pdf("all_env_clines_Vars.pdf", width = 10, height = 7)
par(mfrow=c(2,4))
for(i in c(1:2,4:8))
{hzar.plot.cline(c.fit.env[[i]]$analysis$model.selected, main=colnames(TAB_VARS[,-c(1:2)])[i])}
dev.off()

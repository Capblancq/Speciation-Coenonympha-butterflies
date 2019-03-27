##############################################################################################
##############################################################################################

# Script morphologic clines analyses (Coenonypmha) 
# Capblancq Thibaut
# Date : 31-10-2018

##############################################################################################
##############################################################################################

setwd("./")

library(ggplot2)
library(ggExtra)
library(raster)
library(gridExtra)
library(cowplot)
library(hzar)
library(foreach)
library(doParallel)

##############
#### Data ####

## info
info<-read.table("./infos_ind_VARS.txt", header=T, sep = ",")

## Morphometric traits coming from Script_morphometrics.R
TAB_scaled

## Distance along the hybrid zone coming from Script_Genomic_Cline.R
dist_Vars<-pointDistance(c(997027.1, 6383629.52),info[info$Zone=="VARS",4:5], lonlat=FALSE)
dist_Vars_pop <- aggregate(dist_Vars, by = list(info$Pop[info$Zone=="VARS"]), mean)


##############
#### Vars ####

TAB_VARS <- TAB_scaled[which(TAB_scaled$Pop%in%as.character(dist_Vars_pop$Group.1)),]

#########################
#### Cline with hzar ####

registerDoParallel(cores = 3)
c.fit.traits = foreach (i = colnames(TAB_VARS[,-c(1:2)])) %dopar%
{
  
  c.traits <- list()
  ## Space to hold the observed data
  c.traits$obs <- list()
  ## Space to hold the models to fit
  c.traits$models <- list()
  ## Space to hold the compiled fit requests
  c.traits$fitRs <- list()
  ## Space to hold the output data chains
  c.traits$runs <- list()
  ## Space to hold the analysed data
  c.traits$analysis <- list()
  
  ## Formating data
  c.traits$obs <- hzar.doNormalData1DRaw(hzar.mapSiteDist(dist_Vars_pop$Group.1,dist_Vars_pop$x),TAB_VARS$Pop,TAB_VARS[,i])
  
  ## Looking at the observed data
  # hzar.plot.obsData(c.traits$obs)
  
  ## Creating models 
  c.traits.loadmodel <- function(scaling,tails,id=paste(scaling,tails,sep="."))
  {
    c.traits$models[[id]] <<- hzar.makeCline1DNormal(c.traits$obs, tails)
    if (all(regexpr("fixed",scaling,ignore.case=TRUE) == 1 ))
    {hzar.meta.fix(c.traits$models[[id]])$muL <<- TRUE
    hzar.meta.fix(c.traits$models[[id]])$muR <<- TRUE
    hzar.meta.fix(c.traits$models[[id]])$varL <<- TRUE
    hzar.meta.fix(c.traits$models[[id]])$varR <<- TRUE
    }
    hzar.meta.init(c.traits$models[[id]])$muL <<- c.traits$obs$frame["TOU","mu"]
    #hzar.meta.init(c.traits$models[[id]])$varL <<- c.traits$obs$frame["TOU","var"]
    hzar.meta.init(c.traits$models[[id]])$muR <<- c.traits$obs$frame["VAR_N","mu"]
    #hzar.meta.init(c.traits$models[[id]])$varR <<- c.traits$obs$frame["VAR_N","var"]
  }
  
  c.traits.loadmodel("fixed","none","model.none")
  c.traits.loadmodel("fixed","both","model.both")
  c.traits.loadmodel("fixed","right","model.right")
  c.traits.loadmodel("fixed","left","model.left")
  c.traits.loadmodel("fixed","mirror","model.mirror")
  c.traits.loadmodel("free","none","modelII")
  c.traits.loadmodel("free","both","modelIII")
  
  ## Modifying model parameters
  # Bounds
  c.traits$models <- sapply(c.traits$models,hzar.model.addBoxReq, 0, 9000, simplify=FALSE)
  # reduce the tune setting of modelII from 1.5 to 1.4
  hzar.meta.tune(c.traits$models$modelII)<-1.4
  # reduce the tune setting of modelIII from 1.5 to 1.2
  hzar.meta.tune(c.traits$models$modelIII)<-1.2
  
  ## Compile models to prepare for fitting
  c.traits$fitRs$init <- sapply(c.traits$models, hzar.first.fitRequest.gC, obsData=c.traits$obs, verbose=FALSE, simplify=FALSE)
  
  ## Model chain length and seeds
  chainLength=1e5;                       
  mainSeed=list(A=c(978,544,99,596,528,124), B=c(544,99,596,528,124,978), C=c(99,596,528,124,978,544), D=c(978,99,596,528,124,544), E=c(124,596,978,99,528,544), F=c(99,124,596,978,528,544), G=c(528,99,124,596,978,544))
  lapply(c.traits$fitRs$init, function(x) x$mcmcParam$chainLength <- chainLength)
  lapply(c.traits$fitRs$init, function(x) x$mcmcParam$burnin <- chainLength %/% 10)
  lapply(1:length(c.traits$fitRs$init), function(x) c.traits$fitRs$init[[x]]$mcmcParam$seed[[1]] <- mainSeed[[x]])
  
  ## Run model for an initial chain
  c.traits$runs$init <- list()
  c.traits$runs$init$model.none <- hzar.doFit(c.traits$fitRs$init$model.none)
  c.traits$runs$init$model.both <- hzar.doFit(c.traits$fitRs$init$model.both)
  c.traits$runs$init$model.right <- hzar.doFit(c.traits$fitRs$init$model.right)
  c.traits$runs$init$model.left <- hzar.doFit(c.traits$fitRs$init$model.left)
  c.traits$runs$init$model.mirror <- hzar.doFit(c.traits$fitRs$init$model.mirror)
  c.traits$runs$init$modelII <- hzar.doFit(c.traits$fitRs$init$modelII)
  c.traits$runs$init$modelIII <- hzar.doFit(c.traits$fitRs$init$modelIII)
  
  ## Compile a new set of fi requests using the initial chains
  c.traits$fitRs$chains <- lapply(c.traits$runs$init, hzar.next.fitRequest)
  
  ## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
  c.traits$fitRs$chains <- hzar.multiFitRequest(c.traits$fitRs$chains, each=3, baseSeed=NULL)
  
  ## Just to be thorough, randomize the initial value for each fit
  # center for all models
  lapply(1:length(c.traits$fitRs$chains), function(x) c.traits$fitRs$chains[[x]]$modelParam$init["center"]=runif(21,0,9000)[x])
  # width for all models
  lapply(1:length(c.traits$fitRs$chains), function(x) c.traits$fitRs$chains[[x]]$modelParam$init["width"]=runif(21,0,10000)[x])
  # varH for all models
  lapply(1:length(c.traits$fitRs$chains), function(x) c.traits$fitRs$chains[[x]]$modelParam$init["varH"]=10^runif(21,-1,1)[x])
  
  ## Go ahead and run a chain of 3 runs for every fit request
  c.traits$runs$chains <-  hzar.doChain.multi(c.traits$fitRs$chains, doPar=TRUE, inOrder=FALSE, count=3)
  
  ## Start aggregation of data for analysis
  c.traits$analysis$initDGs <- list()
  
  ## Create a model data group for each model
  c.traits$analysis$initDGs$model.none <- hzar.dataGroup.add(c.traits$runs$init$model.none)
  c.traits$analysis$initDGs$model.both <- hzar.dataGroup.add(c.traits$runs$init$model.both)
  c.traits$analysis$initDGs$model.right <- hzar.dataGroup.add(c.traits$runs$init$model.right)
  c.traits$analysis$initDGs$model.left<- hzar.dataGroup.add(c.traits$runs$init$model.left)
  c.traits$analysis$initDGs$model.mirror <- hzar.dataGroup.add(c.traits$runs$init$model.mirror)
  c.traits$analysis$initDGs$modelII <- hzar.dataGroup.add(c.traits$runs$init$modelII)
  c.traits$analysis$initDGs$modelIII <- hzar.dataGroup.add(c.traits$runs$init$modelIII)
  
  ## Create a hzar.obsDataGroup object from the models
  c.traits$analysis$oDG <- hzar.make.obsDataGroup(c.traits$analysis$initDGs)
  c.traits$analysis$oDG <- hzar.copyModelLabels(c.traits$analysis$initDGs, c.traits$analysis$oDG)
  
  ## Convert all runs to hzar.dataGroup objects
  c.traits$analysis$oDG <- hzar.make.obsDataGroup(lapply(c.traits$runs$chains, hzar.dataGroup.add), c.traits$analysis$oDG)
  
  ## Do model selection based on the AICs scores
  c.traits$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(c.traits$analysis$oDG)
  
  ## Print out the model with the minimum AICc score
  c.traits$analysis$model.name <- rownames(c.traits$analysis$AICcTable)[[which.min(c.traits$analysis$AICcTable$AICc)]]
  
  ## Extract the hzar.dataGroup object for the selected model
  c.traits$analysis$model.selected <- c.traits$analysis$oDG$data.groups[[c.traits$analysis$model.name]]
  
  return(c.traits)
  
}

########################
#### Model selected ####

models <- as.data.frame(lapply(c.fit.traits, function(x) c(unlist(which.min(hzar.AICc.hzar.obsDataGroup(x$analysis$oDG)[,1])), unlist(hzar.AICc.hzar.obsDataGroup(x$analysis$oDG)[which.min(hzar.AICc.hzar.obsDataGroup(x$analysis$oDG)[,1]),]))))
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
for (i in 1:length(c.fit.traits))
{
  c<-density(c.fit.traits[[i]]$analysis$model.selected$data.param$center)
  center<-c$x[which.max(c$y)]
  w<-density(c.fit.traits[[i]]$analysis$model.selected$data.param$width)
  width<-w$x[which.max(w$y)]
  slope<-(abs(c.fit.traits[[i]]$analysis$model.selected$ML.cline$param.all$muR-c.fit.traits[[i]]$analysis$model.selected$ML.cline$param.all$muL)/width)
  par<-c(center, width, slope)
  params_VARS<-rbind(params_VARS, par)
}
row.names(params_VARS)<-colnames(TAB_VARS[,-c(1:2)])
colnames(params_VARS)<-c("center", "width", "Slope")
params_VARS<-as.data.frame(params_VARS)
params_VARS$AIC <- as.numeric(models[2,])
params_VARS$model <- as.numeric(models[1,])
write.table(params_VARS, "params_clines_traits.csv", sep=",")

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
cline.traits.scaled <- list()
cline.traits.scaled[[1]] <- hzar.gen.cline(list(center = c.fit.traits[[1]]$analysis$model.selected$ML.cline$param.all$center, width = c.fit.traits[[1]]$analysis$model.selected$ML.cline$param.all$width, muL=1,muR=0, varH = c.fit.traits[[1]]$analysis$model.selected$ML.cline$param.all$varH, varL = c.fit.traits[[1]]$analysis$model.selected$ML.cline$param.all$varL, varR = c.fit.traits[[1]]$analysis$model.selected$ML.cline$param.all$varR), c.fit.traits[[1]]$analysis$model.selected)
cline.traits.scaled[[2]] <- hzar.gen.cline(list(center = c.fit.traits[[2]]$analysis$model.selected$ML.cline$param.all$center, width = c.fit.traits[[2]]$analysis$model.selected$ML.cline$param.all$width, muL=1,muR=0, varH = c.fit.traits[[2]]$analysis$model.selected$ML.cline$param.all$varH, varL = c.fit.traits[[2]]$analysis$model.selected$ML.cline$param.all$varL, varR = c.fit.traits[[2]]$analysis$model.selected$ML.cline$param.all$varR), c.fit.traits[[2]]$analysis$model.selected)

cols <- brewer.pal(8, "Dark2")[c(1,6)]

pdf("traits_clines_Vars_bounded.pdf", width = 5, height = 5)
plot(1,1, xlim=c(0,8000), ylim=c(0,1), "n", ylab = "Percent characteristic value", xlab = "Distance (m)")
for(i in c(1:2)){hzar.plot.cline2(cline.traits.scaled[[i]], col=cols[i], add=TRUE, lwd=4, lty=1)}
legend(6000,1, legend = colnames(TAB_VARS[,-c(1:2)]), col = cols, lty=1, lwd=4, cex = 0.8)
dev.off()

pdf("all_traits_clines_Vars.pdf", width = 7, height = 4.5)
par(mfrow=c(1,3))
for(i in c(1:3))
{hzar.plot.cline(c.fit.traits[[i]]$analysis$model.selected, main=colnames(TAB_VARS[,-c(1:2)])[i])}
dev.off()


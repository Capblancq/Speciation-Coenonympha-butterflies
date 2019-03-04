#####################################################################
###      Treatment of dadi results for Coenonympha species       ####
#####################################################################

library(ggplot2)

#### Results arcania-gardetta and darwiniana-macromma ####
setwd("~/Documents/ProjetPapillons/These_Hybridation_Coenonympha/Papier2_Hybridation_dans_le_complexe/Genetique/dadi/RESULTS/")

species <- c("arcgar_10", "macdar_10")
table_final <- list()
for(i in species)
{
list <- intersect(list.files(pattern = ".log", recursive=TRUE), list.files(pattern = i, recursive=TRUE))
table <- lapply(list, function(x) read.table(x, sep="\t")[-c(2,3),])
names(table) <- unlist(lapply(lapply(strsplit(list,"/Coenonympha_", "log"), function(x) x[3]), function(x) substr(x, 1,12)))

best_runs <- apply(do.call(cbind, lapply(table, function(x) x[,3])),1,which.min)

table_best <- NULL
for(j in 1:16)
{table_best <- rbind(table_best, table[[best_runs[j]]][j,])}

table_final[[i]] <- table_best
}


## Calculating delta AIC, model score and weigthed AIC
delta_AIC <- lapply(table_final, function(x) x[,3] - x[which.min(x[,3]),3])
mod_score <- lapply(delta_AIC, function(x) (x[which.max(unlist(x))] - x) / x[which.max(unlist(x))])
wAIC <- lapply(delta_AIC, function(x) exp((-x)/2) / sum(exp((-x)/2)))

delta_AIC <- do.call(rbind, delta_AIC)
colnames(delta_AIC) <- table_final[[1]][,1]
mod_score <- do.call(rbind, mod_score)
colnames(mod_score) <- table_final[[1]][,1]
wAIC <- do.call(rbind, wAIC)
colnames(wAIC) <- table_final[[1]][,1]


## Graphical representation
setwd("~/Documents/ProjetPapillons/These_Hybridation_Coenonympha/Papier2_Hybridation_dans_le_complexe/Genetique/dadi/")


pdf("./wAIC_grid.pdf", width = 4, height = 1.5)
grid <- NULL
grid <- expand.grid(x=1:14, y=1:2)
grid$z <- as.vector(t(wAIC))
ggplot(grid, aes(x, y)) + 
         geom_raster(aes(fill=z)) +
         scale_y_continuous(breaks=c(1:2), labels=rownames(wAIC), expand = c(0,0)) + 
         scale_x_continuous(breaks=c(1:14), labels=colnames(wAIC), expand = c(0,0)) +
         xlab("") + ylab("") +
         scale_fill_gradient(low = "yellow", high = "red") +
         theme(axis.text.x = element_text(size=8, hjust=1, angle=60), axis.text.y = element_text(size=8), axis.line = element_blank(), panel.border = element_rect(fill=NA, colour = "black", size=1, linetype="solid"), plot.title = element_text(size=10))
dev.off()

pdf("./delta_AIC_grid.pdf", width = 4, height = 1.5)
grid <- NULL
grid <- expand.grid(x=1:14, y=1:2)
grid$z <- as.vector(t(delta_AIC))
ggplot(grid, aes(x, y)) + 
  geom_raster(aes(fill=z)) +
  scale_y_continuous(breaks=c(1:2), labels=rownames(delta_AIC), expand = c(0,0)) + 
  scale_x_continuous(breaks=c(1:14), labels=colnames(delta_AIC), expand = c(0,0)) +
  xlab("") + ylab("") +
  scale_fill_gradient(low = "yellow", high = "red") +
  theme(axis.text.x = element_text(size=8, hjust=1, angle=60), axis.text.y = element_text(size=8), axis.line = element_blank(), panel.border = element_rect(fill=NA, colour = "black", size=1, linetype="solid"), plot.title = element_text(size=10))
dev.off()

pdf("./mod_score_grid.pdf", width = 4, height = 1.5)
grid <- NULL
grid <- expand.grid(x=1:14, y=1:2)
grid$z <- as.vector(t(mod_score))
ggplot(grid, aes(x, y)) + 
  geom_raster(aes(fill=z)) +
  scale_y_continuous(breaks=c(1:2), labels=rownames(mod_score), expand = c(0,0)) + 
  scale_x_continuous(breaks=c(1:14), labels=colnames(mod_score), expand = c(0,0)) +
  xlab("") + ylab("") +
  scale_fill_gradient(low = "yellow", high = "red") +
  theme(axis.text.x = element_text(size=8, hjust=1, angle=60), axis.text.y = element_text(size=8), axis.line = element_blank(), panel.border = element_rect(fill=NA, colour = "black", size=1, linetype="solid"), plot.title = element_text(size=10))
dev.off()

## Table with retained models (delat AIC < 10) for each couple of species 

table_retained <- NULL
names <- NULL
wAIC_retained <- NULL
deltaAIC_retained <- NULL
for (i in species)
{
  lines <- as.data.frame(table_final[[i]][which(delta_AIC[i,]<10),])
  names <- c(names, rep(i, nrow(lines)))
  deltaAIC_retained <- c(deltaAIC_retained, delta_AIC[i,which(delta_AIC[i,]<10)])
  wAIC_retained <- c(wAIC_retained, wAIC[i,which(delta_AIC[i,]<10)])
  table_retained <- rbind(table_retained, lines)
}
table_retained <- cbind(names, table_retained[,1:3], deltaAIC_retained, wAIC_retained, table_retained[,4:5])
colnames(table_retained) <- c("SPECIES", "MODEL", "MLE", "AIC", "deltaAIC", "wAIC", "theta", "Parameters")

#######################################################################
#######################################################################

### Results triplet

setwd("~/Documents/ProjetPapillons/These_Hybridation_Coenonympha/Papier2_Hybridation_dans_le_complexe/Genetique/dadi/")

dir <- c("RESULTS_HS_fullgeneflow/","RESULTS_HS_nogeneflow/", "RESULTS_HS_geneflow_arcania/", "RESULTS_HS_geneflow_gardetta/", "RESULTS_SGF_arcania/", "RESULTS_SGF_gardetta/")

table_res <- list()
for(DIR in dir)
{
  setwd(DIR)
  species <- c("arcgardar_10", "arcgarmac_10")
  for(i in species)
  {
    list <- intersect(list.files(pattern = ".log", recursive=TRUE), list.files(pattern = i, recursive=TRUE))
    table <- lapply(list, function(x) read.table(x, sep="\t"))
    names(table) <- unlist(lapply(lapply(strsplit(list,"/Coenonympha_", "log"), function(x) x[2]), function(x) substr(x, 1,12)))
    
    best_runs <- apply(do.call(cbind, lapply(table, function(x) x[,2])),1,which.min)
    
    table_best <- NULL
    for(j in 1:16)
    {table_best <- rbind(table_best, table[[best_runs[j]]][j,])}
    
    table_res[[DIR]][[i]] <- table_best
  }
  setwd("~/Documents/ProjetPapillons/These_Hybridation_Coenonympha/Papier2_Hybridation_dans_le_complexe/Genetique/dadi/")
}

## All
TAB_arcgardar <- do.call(rbind,lapply(table_res, function(x) x[[1]]))
TAB_arcgarmac <- do.call(rbind,lapply(table_res, function(x) x[[2]]))

TAB <-list(TAB_arcgardar, TAB_arcgarmac)
names(TAB) <- c('arcgardar', 'arcgarmac')

## Calculating delta AIC, model score and weigthed AIC
delta_AIC <- lapply(TAB, function(x) x[,2] - x[which.min(x[,2]),2])
mod_score <- lapply(delta_AIC, function(x) (x[which.max(unlist(x))] - x) / x[which.max(unlist(x))])
wAIC <- lapply(delta_AIC, function(x) exp((-x)/2) / sum(exp((-x)/2)))

delta_AIC <- do.call(rbind, delta_AIC)
colnames(delta_AIC) <- row.names(TAB[[1]])
mod_score <- do.call(rbind, mod_score)
colnames(mod_score) <- row.names(TAB[[1]])
wAIC <- do.call(rbind, wAIC)
colnames(wAIC) <- row.names(TAB[[1]])


## Figures
pdf("./wAIC_grid_Triplet.pdf", width = 6, height = 4)
grid <- NULL
grid <- expand.grid(x=1:6, y=1:2)
grid$z <- as.vector(t(wAIC))
ggplot(grid, aes(x, y)) + 
  geom_raster(aes(fill=z)) +
  scale_y_continuous(breaks=c(1:2), labels=rownames(wAIC), expand = c(0,0)) + 
  scale_x_continuous(breaks=c(1:6), labels=colnames(wAIC), expand = c(0,0)) +
  xlab("") + ylab("") +
  scale_fill_gradient(low = "yellow", high = "red") +
  theme(axis.text.x = element_text(size=8, hjust=1, angle=60), axis.text.y = element_text(size=8), axis.line = element_blank(), panel.border = element_rect(fill=NA, colour = "black", size=1, linetype="solid"), plot.title = element_text(size=10))
dev.off() 

pdf("./mod_score_grid_Triplet.pdf", width = 6, height = 4)
grid <- NULL
grid <- expand.grid(x=1:6, y=1:2)
grid$z <- as.vector(t(mod_score))
ggplot(grid, aes(x, y)) + 
  geom_raster(aes(fill=z)) +
  scale_y_continuous(breaks=c(1:2), labels=rownames(mod_score), expand = c(0,0)) + 
  scale_x_continuous(breaks=c(1:6), labels=colnames(mod_score), expand = c(0,0)) +
  xlab("") + ylab("") +
  scale_fill_gradient(low = "yellow", high = "red") +
  theme(axis.text.x = element_text(size=8, hjust=1, angle=60), axis.text.y = element_text(size=8), axis.line = element_blank(), panel.border = element_rect(fill=NA, colour = "black", size=1, linetype="solid"), plot.title = element_text(size=10))
dev.off()


## Table scenario retained

species <- c("arcgardar", "arcgarmac")
table_retained_triplet <- NULL
names <- NULL
wAIC_retained <- NULL
deltaAIC_retained <- NULL
for (i in species)
{
  lines <- as.data.frame(TAB[[i]][which(delta_AIC[i,]<10),])
  names <- c(names, rep(i, nrow(lines)))
  deltaAIC_retained <- c(deltaAIC_retained, delta_AIC[i,which(delta_AIC[i,]<10)])
  wAIC_retained <- c(wAIC_retained, wAIC[i,which(delta_AIC[i,]<10)])
  table_retained_triplet <- rbind(table_retained_triplet, lines)
}
table_retained_triplet <- cbind(names, row.names(table_retained_triplet), table_retained_triplet[,1:2], deltaAIC_retained, wAIC_retained, table_retained_triplet[,3:4])
colnames(table_retained_triplet) <- c("SPECIES", "MODEL", "MLE", "AIC", "deltaAIC", "wAIC", "theta", "Parameters")

########################################################################
########################################################################
setwd("~/Documents/ProjetPapillons/These_Hybridation_Coenonympha/Papier2_Hybridation_dans_le_complexe/Genetique/dadi/")

table <- rbind(table_retained, table_retained_triplet)

write.table(x = table, file = "table_retained_dadimodels.csv")



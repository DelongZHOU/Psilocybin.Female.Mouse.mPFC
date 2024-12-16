library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(ggplot2)
library(data.table)
library(Matrix)
library(CellChat)
options(stringsAsFactors = FALSE)
library(unix)
rlimit_as(1e13)

#prepare cell chat db
CellChatDB<-CellChatDB.mouse
CellChatDB.use <- CellChatDB

#load seurat object
combined<-readRDS('/path/to/rds')


#prepare function to initiate cellchat object
initiate_cellchat=function(cells){
Idents(cells)<-cells[['subtype']]
#prepare count matrix
data.input <- GetAssayData(cells, assay = "RNA", slot = "data") #extract normalized data matrix
labels <-factor(Idents(cells),levels=lvl) #extract identity/label
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group") #create cellchat object
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels

#start calculating communication network
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)
return(cellchat)
}

#subset seurat object by condition
sal<-subset(combined,subset = sex=='F' & condition=='SAL')
psy<-subset(combined,subset = sex=='F' & condition=='PSY')

#make cellchat objects
cellchat_sal<-initiate_cellchat(sal)
cellchat_psy<-initiate_cellchat(psy)

#merge cellchat
object.list<-c(cellchat_sal,cellchat_psy)
names(object.list)<-c('Sal','Psy')
cellchat <- mergeCellChat(object.list, add.names = c('Sal','Psy'))

#save count and table
for (i in 1:length(object.list)) {
file.name<-'/path/to/csv'
write.csv(object.list[[i]]@net$count,file.name)
}

#circle plot for difference in interaction number
#reverse color as the diff is backwards
png('/path/to/png',width = 2400, height = 2400, units = "px", res=600)
netVisual_diffInteraction(cellchat,comparison = c(1, 2),top = 0.10, weight.scale = T,arrow.size = 1,edge.width.max = 4, color.edge = c('#E6308A','#808080'),color.use=c('#c6dbef','#9ecae1','#3182bd','#08519c','#33a02c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'))
dev.off()
pdf('/path/to/pdf')
netVisual_diffInteraction(cellchat,comparison = c(1, 2),top = 0.10, weight.scale = T,arrow.size = 1, edge.width.max = 4, color.edge = c('#E6308A','#808080'),color.use=c('#c6dbef','#9ecae1','#3182bd','#08519c','#33a02c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'))
dev.off()

#circle plot for difference in interaction strength
#reverse color as the diff is backwards
png('/path/to/png',width = 2400, height = 2400, units = "px", res=600)
netVisual_diffInteraction(cellchat,comparison = c(1, 2),top = 0.10, measure = "weight", weight.scale = T,arrow.size = 1,edge.width.max = 4, color.edge = c('#E6308A','#808080'),color.use=c('#c6dbef','#9ecae1','#3182bd','#08519c','#33a02c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'))
dev.off()
pdf('/path/to/pdf')
netVisual_diffInteraction(cellchat,comparison = c(1, 2),top = 0.10, measure = "weight", weight.scale = T,arrow.size = 1, edge.width.max = 4, color.edge = c('#E6308A','#808080'),color.use=c('#c6dbef','#9ecae1','#3182bd','#08519c','#33a02c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'))
dev.off()

library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(ggplot2)
library(data.table)
library(Matrix)
library(harmony)
library(tidyr)
library(unix)
rlimit_as(1e13)

#read individual seurat files
samples=list.files('/path/to/seurat')
seurat_array = lapply(samples,readRDS)
#merge seurat objects
seurat_tail = tail(seurat_array,-1)
combined = merge(seurat_array[[1]], seurat_tail)

#UMAP
combined <- subset(combined, nFeature_RNA > 800)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst")
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)
combined <- RunUMAP(combined, reduction =  "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction =  "pca", dims = 1:30)
saveRDS(combined, '/path/to/rds')

#QC plot
p1 <- DimPlot(combined, reduction = "umap",label = TRUE)+ NoLegend()
ggsave(file="psy.F.combined.dim.umap.png",width=7,height=7)
p1 <- FeaturePlot(combined,features = 'nFeature_RNA', reduction = "umap",max.cutoff = 2000)
ggsave(file="psy.F.combined.dim.umap.nFeature_RNA.png",width=7,height=7)
p1 <- FeaturePlot(combined,features = 'nCount_RNA', reduction = "umap",max.cutoff = 20000)
ggsave(file="psy.F.combined.dim.umap.nCount_RNA.png",width=7,height=7)


#cell type marker plot
cell.type.markers=c('Gja1','Flt1','C1qa','Bmp4','Aspa','Pdgfra','Snap25','Slc17a7','Gad2','Sdk1','Sox11')

#try different resolutions for clustering
for (res in c(0.1,0.2,0.5,1,2)) {
combined <- FindClusters(combined, resolution = res)
cluster_labels<-table(combined[[c('seurat_clusters','label.from.pfc','label.from.allen','label')]])
file.title=paste0('psy.F.combined.umap.res=',res,'.csv')
write.csv(cluster_labels,file.title)
p1 <- DimPlot(combined, reduction = "umap", label = TRUE,shuffle=TRUE)+ NoLegend()
file.title=paste0('psy.F.combined.umap.res=',res,'.png')
ggsave(file=file.title,width=7,height=7)
p1=VlnPlot(combined,'nFeature_RNA',pt.size = 0) + NoLegend()
file.title=paste0('psy.F.combined.nfeature.vln.res=',res,'.png')
ggsave(file.title,height=7,width=7)
file.title=paste0('psy.F.combined.nCount_RNA.vln.res=',res,'.png')
p1=VlnPlot(combined,'nCount_RNA',pt.size = 0) + NoLegend()
ggsave(file.title,height=7,width=7)
file.title=paste0('psy.F.combined.percent.mt.vln.res=',res,'.png')
p1=VlnPlot(combined,'percent.mt',pt.size = 0) + NoLegend()
ggsave(file.title,height=7,width=7)
file.title=paste0('psy.F.combined.cell.type.markers.vln.res=',res,'.png')
p1=VlnPlot(combined,cell.type.markers,pt.size = 0, stack=TRUE) + NoLegend()
ggsave(file.title,height=7,width=7)
}

#After inspection, resolution=0.5 is selected
combined <- FindClusters(combined, resolution = 0.5)

#further QC and subclustering for neurons
#select the clusters surviving the QC
gluta_list<-c(5,7,17,18)
gluta_name<-c('L2-6 IT','L6 CT','L5-6 NP','L5 PT')
names(gluta_name)<-gluta_list
gluta<-subset(combined,idents=gluta_list)
gluta<-RenameIdents(gluta,gluta_name)
gluta[['subtype']]<-Idents(gluta)
gluta[['cell.type']]<-'Gluta'
gluta<-subset(gluta, subset = nFeature_RNA > 800)

cell.type.markers=c('Gja1','Flt1','C1qa','Bmp4','Aspa','Pdgfra','Snap25','Slc17a7','Gad2','Sdk1')
gluta.marker<-c("Cux2","Calb1","Etv1","Syt6")
gluta <- NormalizeData(gluta)
gluta <- FindVariableFeatures(gluta, selection.method = "vst")
gluta <- ScaleData(gluta)
gluta <- RunPCA(gluta, npcs = 30)
gluta <- RunUMAP(gluta, dims = 1:30)

p1 <- FeaturePlot(gluta,features = 'nFeature_RNA', reduction = "umap",max.cutoff = 4000)
ggsave(file="psy.gluta.umap.nFeature_RNA.png",width=7,height=7)
p1 <- FeaturePlot(gluta,features = 'nCount_RNA', reduction = "umap",max.cutoff = 20000)
ggsave(file="psy.gluta.umap.nCount_RNA.png",width=7,height=7)
p1 <- DimPlot(gluta, reduction = "umap", label = TRUE, ,shuffle=TRUE)+ NoLegend()
ggsave(file="psy.gluta.umap.label.png",width=7,height=7)

gluta <- FindNeighbors(gluta, reduction =  "pca", dims = 1:30)
for (res in c(0.1,0.2,0.5)) {
gluta <- FindClusters(gluta, resolution = res)
p1 <- DimPlot(gluta, reduction = "umap", label = TRUE, ,shuffle=TRUE)+ NoLegend()
file.title=paste0('psy.gluta.umap.res=',res,'.png')
ggsave(file=file.title,width=7,height=7)
p1=VlnPlot(gluta,cell.type.markers,pt.size = 0, stack=TRUE) + NoLegend()
file.title=paste0('psy.gluta.celltype.markers.vln.res=',res,'.png')
ggsave(file.title,height=7,width=7)
p1=VlnPlot(gluta,gluta.marker,pt.size = 0, stack=TRUE) + NoLegend()
file.title=paste0('psy.gluta.subtype.markers.vln.res=',res,'.png')
ggsave(file.title,height=7,width=7)
p1=VlnPlot(gluta,'nFeature_RNA',pt.size = 0) + NoLegend()
file.title=paste0('psy.gluta.nfeature.vln.res=',res,'.png')
ggsave(file.title,height=7,width=7)
file.title=paste0('psy.gluta.nCount_RNA.vln.res=',res,'.png')
p1=VlnPlot(gluta,'nCount_RNA',pt.size = 0) + NoLegend()
ggsave(file.title,height=7,width=7)
}
#After inspection, resolution=0.1 is selected
gluta <- FindClusters(gluta, resolution = 0.1)
gluta_list<-c(0,1,2,3)
gluta_name<-c('L2-6 IT','L6 CT','L5-6 NP','L5 PT')
names(gluta_name)<-gluta_list
gluta<-subset(combined,idents=gluta_list)
gluta<-RenameIdents(gluta,gluta_name)
gluta[['subtype']]<-Idents(gluta)


#GABA cells are too few to subcluster
gaba<-subset(combined,idents=c(13))
gaba[['subtype']]<-'GABA'
gaba[['cell.type']]<-'GABA'
gaba<-subset(gaba, subset = nFeature_RNA > 800)

#Glia
glia_list<-c(1,2,8,4,19,6,10,0,3)
glia_name<-c('Astro','Astro','Astro','Endo','Endo','Micro','NFO','Oligo','OPC')
names(glia_name)<-glia_list
glia<-subset(combined,idents=glia_list)
glia<-RenameIdents(glia,gluta_name)
glia[['subtype']]<-Idents(glia)
glia[['cell.type']]<-Idents(glia)

combined = merge(glia,c(gluta,gaba))
saveRDS(combined,'path/to/rds')


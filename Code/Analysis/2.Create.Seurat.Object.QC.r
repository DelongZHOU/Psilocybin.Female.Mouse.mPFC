library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(ggplot2)
library(data.table)
library(Matrix)
library(Libra)
library(harmony)
library(CellChat)
options(stringsAsFactors = FALSE)
library(unix)
rlimit_as(1e13)

samples=list.files('/path/to/cellbender')

create_seurat = function(name){
#determine condition and sex from file name
condition=substr(name,1,3)
sex=substr(name,5,5)

count.path=paste0('/path/to/cellbender',name)
count=Read10X_h5(count.path, use.names = TRUE)
seurat_object=CreateSeuratObject(count,min.cells = 3, min.features = 200)
name<-substr(name,1,6)
seurat_object[["sample"]]=name
seurat_object[['condition']]=condition
seurat_object[['sex']]=sex

#pre-process
#nFeature
nFeature<-seurat_object$nFeature_RNA
average<-mean(nFeature)
stdev<-sd(nFeature)
nFeature_min<-max(0,average-2*stdev)
nFeature_max<-average+2*stdev
#nCount
nCount<-seurat_object$nCount_RNA
average<-mean(nCount)
stdev<-sd(nCount)
nCount_min<-max(0,average-2*stdev)
nCount_max<-average+2*stdev
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
seurat_object<-subset(seurat_object, nCount_RNA > nCount_min & 
									 nCount_RNA < nCount_max & 
									 nFeature_RNA > nFeature_min & 
									 nFeature_RNA < nFeature_max & 
									 nFeature_RNA > 200 & 
									 percent.mt < 10
									 )

file.title<-paste0('/path/to/seurat/')
saveRDS(seurat_object,file.title)
}

lapply(samples,create_seurat)
library(dplyr)
library(Seurat)
library(ggplot2)

#load data and order cell type
combined<-readRDS('path/to/rds')
Idents(combined)<-combined$subtype

cell.type.markers=c('Snap25','Slc17a7','Cux2','Bcl6','Tshz2','Syt6','Gad2','Gja1','Flt1','C1qa','Bmp4','Aspa','Pdgfra')
cell.type.marker.color=c('#7fcdbb','#6baed6','#c6dbef','#9ecae1','#3182bd','#08519c','#33a02c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
cell.type.color=c('#c6dbef','#9ecae1','#3182bd','#08519c','#33a02c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')

p1 <- DimPlot(combined, reduction = 'umap',label = TRUE,cols=cell.type.color)+ NoLegend()
ggsave('Fig.1.H.UMAP.png',height=4,width=4,dpi=600)
ggsave('Fig.1.H.UMAP.pdf',height=4,width=4,dpi=600)

p1=VlnPlot(combined,cell.type.markers,pt.size = 0, stack=TRUE,cols=cell.type.marker.color) + NoLegend()
ggsave('Fig.1.I.celltype.markers.png',height=4,width=4,dpi=600)
ggsave('Fig.1.I.celltype.markers.pdf',height=4,width=4)

library(dplyr)
library(Seurat)
library(ggplot2)
FontSize(12)

#load data and order cell type
combined<-readRDS('path/to/rds')
#vlnplt qc
Idents(combined)<-combined$sample_id
p1 <- VlnPlot(combined,features = 'percent.mt',pt.size = 0,cols=c('#808080','#808080','#808080','#808080','#E6308A','#E6308A','#E6308A','#E6308A')) + NoLegend() + geom_boxplot(width = 0.2, outlier.shape = NA) + ggtitle('Mito Percent')
ggsave("Fig.S1.A.sample.percent.mt.png",width=3,height=3)
ggsave("Fig.S1.A.sample.percent.mt.pdf",width=3,height=3)
p1 <- VlnPlot(combined,features = 'nCount_RNA',pt.size = 0,cols=c('#808080','#808080','#808080','#808080','#E6308A','#E6308A','#E6308A','#E6308A')) + NoLegend() + geom_boxplot(width = 0.2, outlier.shape = NA) + ggtitle('Unique UMIs')
ggsave("Fig.S1.B.sample.total.umi.png",width=3,height=3)
ggsave("Fig.S1.B.sample.total.umi.pdf",width=3,height=3)
p1 <- VlnPlot(combined,features = 'nFeature_RNA',pt.size = 0,cols=c('#808080','#808080','#808080','#808080','#E6308A','#E6308A','#E6308A','#E6308A')) + NoLegend() + geom_boxplot(width = 0.2, outlier.shape = NA) + ggtitle('Unique Genes')
ggsave("Fig.S1.C.sample.unique.genes.png",width=3,height=3)
ggsave("Fig.S1.C.sample.unique.genes.pdf",width=3,height=3)
Idents(combined)<-combined$cell.type
p1 <- VlnPlot(combined,features = 'nFeature_RNA',pt.size = 0,cols=c('#6baed6','#33a02c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')) + NoLegend() + geom_boxplot(width = 0.2, outlier.shape = NA) + ggtitle('Unique Genes')
ggsave("Fig.S1.D.celltype.unique.genes.png",width=3,height=3)
ggsave("Fig.S1.D.celltype.unique.genes.pdf",width=3,height=3)

#order cell types
Idents(combined)<-combined$subtype
cell.type.color=c('#c6dbef','#9ecae1','#3182bd','#08519c','#33a02c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
p1 <- DimPlot(combined, reduction = 'umap',label = FALSE,cols=cell.type.color, split.by='sample_id', ncol=4)+ NoLegend()
ggsave('Fig.S1.F.UMAP.sample.png',height=3,width=4,dpi=600)
ggsave('Fig.S1.F.UMAP.sample.pdf',height=3,width=4,dpi=600)

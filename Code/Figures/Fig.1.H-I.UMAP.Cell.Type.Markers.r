library(dplyr)
library(Seurat)
library(ggplot2)

#load data and order cell type
combined<-readRDS('path/to/rds')
Idents(combined)<-combined$subtype

cell.type.markers=c('Snap25','Slc17a7','Cux2','Bcl6','Tshz2','Syt6','Gad2','Gja1','Flt1','C1qa','Bmp4','Aspa','Pdgfra')
cell.type.marker.color=c('#7fcdbb','#6baed6','#c6dbef','#9ecae1','#3182bd','#08519c','#33a02c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
cell.type.color=c('#c6dbef','#9ecae1','#3182bd','#08519c','#33a02c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')

p1 <- DimPlot(combined, reduction = 'umap',label = False,cols=cell.type.color)+ NoLegend() + labs(title = NULL, x = NULL, y = NULL) + NoAxes()
ggsave('umap.pdf',height=4,width=4,dpi=600)
ggsave('umap.png',height=4,width=4,dpi=600)

p1=VlnPlot(combined,cell.type.markers,pt.size = 90, stack=TRUE,cols=cell.type.marker.color) + NoLegend() +
  theme( 
    strip.text = element_text(face = "plain"), # Remove bold from gene labels
    axis.title.x = element_blank(),            # Remove x-axis label
    axis.ticks.x = element_blank(),            # Remove x-axis ticks
    axis.text.x = element_blank(),             # Remove x-axis range values
    axis.line = element_line(linewidth = 0.5), # Set axis line thickness
    axis.ticks = element_line(linewidth = 0.5) # Set axis line thickness
  )
ggsave('markers.png',height=4,width=4,dpi=600)
ggsave('markers.pdf',height=4,width=4,dpi=600)

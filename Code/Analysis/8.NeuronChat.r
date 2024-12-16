library(Seurat)
library(NeuronChat)
library(CellChat)
library(ggplot2)

#set color
color.preset=c('#808080','#E6308A') #grey=saline pink=psy

#create NeuronChat objects
neuron<-readRDS('/path/to/rds')
subtypes<-sort(unique(neuron$subtype))
neuron$subtype<-factor(neuron$subtype,levels=subtypes)
psy<-subset(neuron,subset=condition=='PSY')
sal<-subset(neuron,subset=condition=='SAL')
create_nc<-function(y){
expr=GetAssayData(y, assay = "RNA", slot = "data")
metadata=y[['subtype']]
x<-createNeuronChat(expr,meta=metadata,DB='mouse',group.by=metadata$subtype)
x <- run_NeuronChat(x,M=100)
return(x)
}
nc_list<-lapply(c(sal,psy),create_nc)
names(nc_list)<-c('SAL','PSI')

#aggregate interaction
net_aggregated_sal_count <- net_aggregation(nc_list[[1]]@net,method='count')
net_aggregated_psy_count <- net_aggregation(nc_list[[2]]@net,method='count')
net_aggregated_sal_weight <- net_aggregation(nc_list[[1]]@net,method='weight')
net_aggregated_psy_weight <- net_aggregation(nc_list[[2]]@net,method='weight')

write.csv(net_aggregated_sal_count,'path/to/csv')
write.csv(net_aggregated_psy_count,'path/to/csv')
write.csv(net_aggregated_sal_weight,'path/to/csv')
write.csv(net_aggregated_psy_weight,'path/to/csv')

#Circle diff plot
#work around using CellChat function

circle_plot<-function(sal,psy,top){
	color.edge = rev(color.preset)
	net <- psy-sal
	net[abs(net) < stats::quantile(abs(net), probs = 1-top, na.rm= T)] <- 0
	g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
	color.use = c('#33a02c','#c6dbef','#9ecae1','#3182bd','#08519c')
	edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
	coords<-layout_(g,in_circle())
	if(nrow(coords)!=1){
	coords_scale=scale(coords)
	}else{
	coords_scale<-coords
	}

	vertex.weight <- 20
	igraph::V(g)$color<-color.use[igraph::V(g)]
	igraph::V(g)$size<-vertex.weight+10
	igraph::V(g)$color<-color.use[igraph::V(g)]
	igraph::V(g)$frame.color <- color.use[igraph::V(g)]
	igraph::V(g)$label.color <- 'black'
	igraph::V(g)$label.cex<-1
	igraph::E(g)$arrow.width<-2
	igraph::E(g)$arrow.size<-0.5
	igraph::E(g)$label.color<-'black'
	igraph::E(g)$label.cex<-0.8
	igraph::E(g)$label<-''
	igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1],color.edge[2])
	igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, 0.6)
	igraph::E(g)$weight <- abs(igraph::E(g)$weight)
	edge.weight.max <- max(igraph::E(g)$weight)
	igraph::E(g)$width<- 0.8+igraph::E(g)$weight/edge.weight.max*0.8
	igraph::E(g)$loop.angle <- 0
	radian.rescale <- function(x, start=0, direction=1) {
		c.rotate <- function(x) (x + start) %% (2 * pi) * direction
		c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
	}
	label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
	label.dist <- vertex.weight/max(vertex.weight)+2
	plot(g,edge.curved=0.2,vertex.shape='circle',layout=coords_scale,margin=1, vertex.label.dist=label.dist,
		   vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
	gg <- recordPlot()
	return(gg)
}

pdf('path/to/pdf')
circle_plot(net_aggregated_sal_count,net_aggregated_psy_count,top=0.1)
dev.off()
png('path/to/png')
circle_plot(net_aggregated_sal_count,net_aggregated_psy_count,top=0.1)
dev.off()
pdf('path/to/pdf')
circle_plot(net_aggregated_sal_weight,net_aggregated_psy_weight,top=0.1)
dev.off()
png('path/to/png')
circle_plot(net_aggregated_sal_weight,net_aggregated_psy_weight,top=0.1)
dev.off()

library(tidyverse)
library(DOSE)
library(clusterProfiler)

#GO flavor
flavors=c('BP','CC','MF')
genes=read.csv('path/to/csv')
cells=unique(genes$cell_type)
for (cell in cells) {
print(cell)
target=genes[genes$cell_type==cell,]$avg_logFC
target_names=genes[genes$cell_type==cell,]$Entrez
names(target)=as.character(target_names)
target<-sort(target,decreasing=TRUE)
#GO
for (flavor in flavors) {
print(flavor)
ego<-gseGO(
target,
OrgDb='org.Mm.eg.db',
ont = flavor,
pvalueCutoff = 0.05,
pAdjustMethod = "BH",
verbose = FALSE
)
write.csv(ego,'path/to/csv')
}
}

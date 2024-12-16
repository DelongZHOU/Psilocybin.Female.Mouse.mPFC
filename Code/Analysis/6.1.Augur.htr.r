library(Seurat)
library(ggplot2)
library(Augur)
library(unix)
rlimit_as(1e13)

#get Htr id from command line
htr_id<- commandArgs(trailingOnly=TRUE)[1]

#these parameters are chosen after a parameter sweeping with vquant={0.5,0.7,0.9}, ntree={50,100,150}, mtry={10,15,20}
#lower vquant in general produces lowest scores and takes more computational time

vquant=0.9
ntree=150
mtry=20

neuron<-readRDS('path/to/rds')

#reduce to htr positive cells
expr <- FetchData(object = neuron, vars = htr_id)
htr_pos=neuron[, which(x = expr > 0)]

htr_pos$cell_type<-factor(htr_pos$subtype)
htr_pos$label<-factor(htr_pos$condition)
expr=GetAssayData(htr_pos,slot='count')
meta = htr_pos@meta.data
augur<-calculate_auc(input = expr, meta = meta, n_threads = 40,var_quantile = vquant,subsample_size = 50,n_subsamples=100,rf_params = list(trees = ntree, mtry = mtry, min_n = NULL, importance = "accuracy"))
write.csv(augur$AUC,'path/to/csv')
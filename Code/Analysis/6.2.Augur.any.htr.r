library(Seurat)
library(ggplot2)
library(Augur)
library(unix)
rlimit_as(1e13)

#these parameters are chosen after a parameter sweeping with vquant={0.5,0.7,0.9}, ntree={50,100,150}, mtry={10,15,20}
#lower vquant in general produces lowest scores and takes more computational time

vquant=0.9
ntree=150
mtry=20

neuron<-readRDS('path/to/rds')

#list of expressed Htr are obtained through examining the expression matrix
htr_pos=subset(neuron,	Htr1a>0 | 
						Htr1b>0 | 
						Htr1d>0 | 
						Htr1f>0 | 
						Htr2a>0 | 
						Htr2b>0 | 
						Htr2c>0 | 
						Htr3a>0 | 
						Htr4>0 | 
						Htr5a>0 | 
						Htr5b>0 | 
						Htr6>0 | 
						Htr7>0
				)
				
Idents(htr_pos)<-factor(htr_pos$subtype)
Idents(htr_neg)<-factor(htr_neg$subtype)

htr_pos$cell_type<-factor(Idents(htr_pos))
htr_pos$label<-factor(htr_pos$condition)
expr=GetAssayData(htr_pos,slot='count')
meta = htr_pos@meta.data
augur<-calculate_auc(input = expr, meta = meta, n_threads = 40,var_quantile = vquant,subsample_size = 50,n_subsamples=100,rf_params = list(trees = ntree, mtry = mtry, min_n = NULL, importance = "accuracy"))
write.csv(augur$AUC,'path/to/csv')

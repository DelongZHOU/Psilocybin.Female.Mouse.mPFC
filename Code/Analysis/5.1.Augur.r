library(Seurat)
library(ggplot2)
library(Augur)
library(unix)
rlimit_as(1e13)

#these parameters are chosen after a parameter sweep with vquant={0.5,0.7,0.9}, ntree={50,100,150}, mtry={10,15,20}
#lower vquant in general produces lowest scores and takes more computational time

vquant=0.9
ntree=150
mtry=20

combined<-readRDS('path/to/rds')
#unfactorized labels used to cause errors
#should be fixed in the latest version
combined$label<-factor(combined$condition)
combined$cell_type<-factor(combined$subtype)

#taking raw counts as suggested by https://github.com/neurorestore/Augur/issues/21
#the default behavior is to use the normalized counts
expr=GetAssayData(combined,slot='count')
meta = combined@meta.data

augur<-calculate_auc(input = expr, meta = meta, n_threads = 40,var_quantile = vquant,subsample_size = 50,n_subsamples=100,rf_params = list(trees = ntree, mtry = ntry, min_n = NULL, importance = "accuracy"))
#the AUC is the model score
write.csv(augur$AUC,'path/to/csv')
#the results contain sampling results
write.csv(augur$results,'path/to/csv')
#use permute to shuffle cells and emperically produce p value
augur<-calculate_auc(input = expr, meta = meta, n_threads = 40,augur_mode = "permute",var_quantile = vquant,subsample_size = 50,n_subsamples=1000,rf_params = list(trees = ntree, mtry = ntry, min_n = NULL, importance = "accuracy"))
write.csv(augur$results,'path/to/csv')

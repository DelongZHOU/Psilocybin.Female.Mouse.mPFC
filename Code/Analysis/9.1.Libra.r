library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(ggplot2)
library(data.table)
library(Matrix)
library(Libra)
library(magrittr)
library(unix)
rlimit_as(1e13)

#load cell type / subtype info from command line
cell_type<- commandArgs(trailingOnly=TRUE)[1]


combined<-readRDS('/path/to/rds')
#subset the data to the cell type in question
combined<-subset(combined,subset=subtype==cell_type)

combined$condition<-factor(combined$condition,levels=c('SAL','PSY'))
DE_f = run_de(combined, de_family = 'mixedmodel', de_method = 'negbinom_offset', de_type = 'Wald', replicate_col = "sample", cell_type_col = "subtype", label_col = "condition",min_reps = 3, n_threads = 40)

write.csv(DE_f,'path/to/csv')
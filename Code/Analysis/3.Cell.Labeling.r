library(Seurat)
library(ggplot2)
library(scibet)
library(tidyverse)
library(viridis)
library(SingleR)
library(SingleCellExperiment)
library(scuttle)
library(singleCellNet)


#scibet
assign.scibet<-function(test,ref){
#extract query exp matrix
query.expr<-GetAssayData(test,slot='data')
query.expr <- as_tibble(t(as.matrix(query.expr)),rownames = NA)
#extract ref exp and label
ref.expr<-GetAssayData(ref,slot='data')
ref.expr<- as_tibble(t(as.matrix(ref.expr)),rownames = NA)
ref.expr<-cbind(ref.expr,label=ref[['label']])

#scibet prediction
scibet.label <- SciBet(train = ref.expr, test = query.expr)
return(scibet.label)
}

#singlecellnet
assign.singlecellnet<-function(test,ref){
#extract test metadata, exp matrix and gene list
query.metadata<-test[[]]
query.exp<-GetAssayData(test,slot='count')
query.genes = rownames(query.exp)
#extract ref metadata, exp matrix and reduce to genes shared between query and ref
ref.metadata<-ref[[]]
ref.exp<-GetAssayData(ref,slot='count')
#filter reference to genes expressed in both set
commonGenes = intersect(rownames(ref.exp), query.genes)
ref.exp.filter = ref.exp[commonGenes,]
#train model from ref
class_info<-scn_train(stTrain = ref.metadata, expTrain = ref.exp.filter, nTopGenes = 10, nRand = 70, nTrees = 1000, nTopGenePairs = 25, dLevel = 'label')
#apply to test
query.scn<-scn_predict(class_info[['cnProc']], query.exp, nrand=0)
#extract predicted label
pred_labels <- assign_cate(classRes = query.scn, sampTab = query.metadata, cThresh = 0.5)
return(pred_labels$category)
}

#singleR
assign.singleR<-function(test,ref){
test <- as.SingleCellExperiment(test)
ref <- as.SingleCellExperiment(ref)
test <- logNormCounts(test) 
ref <- logNormCounts(ref)
preds <- SingleR(test=test, ref=ref, labels=ref$label, de.method="wilcox")
return(preds$labels)
}

#take the majority of cell label from 3 tools
vote<-function(x){
if ((x[1]==x[2]) | (x[1]==x[3])){
return(x[1])
} else {
	if (x[2]==x[3]){return(x[2])} else {return(NA)}
}
}


#load query
name = commandArgs(trailingOnly=TRUE)[1]
file.title<-paste0("/path/to/rds/",name,)
query<-readRDS(file.title)

#load pfc data
pfc<-readRDS('/path/to/Bhattacherjee.rds')
#predict based on addiction set
sal.scibet<-assign.scibet(query,pfc)
sal.singleR<-assign.singleR(query,pfc)
sal.scn<-assign.singlecellnet(query,pfc)

query.labels<-cbind(scibet=sal.scibet,scn=sal.scn,singleR=sal.singleR)
query.vote<-apply(query.labels,1,vote)
query[['label.from.pfc']]<-query.vote


#load allen data
allen<-readRDS('/path/to/yao.rds')
#predict based on allen
sal.scibet<-assign.scibet(query,allen)
sal.singleR<-assign.singleR(query,allen)
sal.scn<-assign.singlecellnet(query,allen)

query.labels<-cbind(scibet=sal.scibet,scn=sal.scn,singleR=sal.singleR)
query.vote<-apply(query.labels,1,vote)
query[['label.from.allen']]<-query.vote


#parse labels
two.labels<-query[[c("label.from.pfc","label.from.allen")]]
sum.label<-function(x){
if (x[1]=='Excitatory'){
	if (x[2] %in% c('Car3','CR','L2_3_IT_CTX','L2_3_IT_ENTl','L2_3_IT_RHP','L2_IT_ENTl','L4_5_IT_CTX','L5_6_IT_TPE-ENT','L5_6_NP_CTX','L5_IT_CTX','L5_PT_CTX','L6_CT_CTX','L6_IT_CTX','L6b_CTX','NP_SUB')) { return(x[2]) } else {return(NA)}
} else if (x[1]=='Inhibitory') {
	if (x[2] %in% c('Lamp5','Meis2','Pvalb','Sncg','Sst','Sst_Chodl','Vip')) { return(x[2]) } else {return(NA)}
} else {
return(x[1])
}
}
labels.use<-drop_na(two.labels)
labels.use['pfc.allen.label']<-apply(labels.use,1,sum.label)
labels.use<-drop_na(labels.use)
cells.use<-rownames(labels.use)
query<-query[,cells.use]
query[['label']]<-query.use['pfc.allen.label']

#save result
query.label.title<-'path/to/rds'
saveRDS(query,query.label.title)
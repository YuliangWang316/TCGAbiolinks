library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(Seurat)
library(dplyr)
library(patchwork)
setwd("E:/GDC/")
clinical<-GDCquery_clinic(project = "CPTAC-3", type = "clinical")
query<-GDCquery(project = "CPTAC-3", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "CellRanger - 10x Raw Counts")
GDCdownload(query)
expdat<-GDCprepare(query = query)





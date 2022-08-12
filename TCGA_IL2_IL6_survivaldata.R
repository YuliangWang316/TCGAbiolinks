library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(survival)
library(survminer)
request_cancer=c("ACC","CHOL","BLCA","COAD","READ","BRCA","LGG","GBM","PCPG","CESC","ESCA","STAD","UVM","HNSC","KIRC","KIRP","KICH","LIHC","LUAD","LUSC","DLBC","LAML","OV","PAAD","MESO","PRAD","SKCM","SARC","TGCT","THYM","THCA","UCEC","UCS")#,"CHOL","BLCA","COAD","READ","BRCA","LGG","GBM","PCPG","CESC","ESCA","STAD","UVM","HNSC","KIRC","KIRP","KICH","LIHC","LUAD","LUSC","DLBC","LAML","OV","PAAD","MESO","PRAD","SKCM","SARC","TGCT","THYM","THCA","UCEC","UCS"
for (i in request_cancer) {
  cancer_type=paste("TCGA",i,sep="-")
  print(cancer_type)
  #下载临床数据
  clinical <- GDCquery_clinic(project = cancer_type, type = "clinical")
  #write.csv(clinical,file = paste(cancer_type,"clinical.csv",sep = "-"))
  
  #下载rna-seq的counts数据
  query <- GDCquery(project = cancer_type, 
                    data.category = "Transcriptome Profiling", 
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts")
  
  GDCdownload(query, method = "api", files.per.chunk = 100)
  expdat <- GDCprepare(query = query)
  count_matrix=assay(expdat)
  
  samplesTP <- TCGAquery_SampleTypes(colnames(count_matrix), typesample = c("TP"))
  IL2 <- count_matrix[c("ENSG00000109471"),samplesTP]
  
  names(IL2) <- sapply(strsplit(names(IL2),'-'),function(x) paste0(x[1:3],collapse="-"))
  IL2<-t(IL2)
  clinical$"IL2" <- IL2[match(clinical$submitter_id,colnames(IL2))]
  IL6 <- count_matrix[c("ENSG00000136244"),samplesTP]
  #IL6 <- as.data.frame(t(IL6))
  names(IL6) <- sapply(strsplit(names(IL6),'-'),function(x) paste0(x[1:3],collapse="-"))
 
  IL6<-t(IL6)
  clinical$"IL6" <- IL6[match(clinical$submitter_id,colnames(IL6))]
  df<-subset(clinical,select =c(submitter_id,vital_status,days_to_death,days_to_last_follow_up,IL2,IL6))
  
  df <- df[!is.na(df$IL2),]
  df <- df[!is.na(df$IL6),]
  df$IL2_status <- ''
  df[df$IL2 > median(df$IL2),]$IL2_status <- "H"
  df[df$IL2 <= median(df$IL2),]$IL2_status <- "L"
  df$IL6_status <- ''
  df[df$IL6 > median(df$IL6),]$IL6_status <- "H"
  df[df$IL6 <= median(df$IL6),]$IL6_status <- "L"
  
  df_IL2_hi<-filter(df,IL2_status == "H")
  df_IL2_hi$IL6_status <- ''
  df_IL2_hi[df_IL2_hi$IL6 > median(df_IL2_hi$IL6),]$IL6_status <- "H"
  df_IL2_hi[df_IL2_hi$IL6 <= median(df_IL2_hi$IL6),]$IL6_status <- "L"
  
  df_IL2_lo<-filter(df,IL2_status == "L")
  df_IL2_lo$IL6_status <- ''
  df_IL2_lo[df_IL2_lo$IL6 > median(df_IL2_lo$IL6),]$IL6_status <- "H"
  df_IL2_lo[df_IL2_lo$IL6 <= median(df_IL2_lo$IL6),]$IL6_status <- "L"
  
  
  df_IL6_hi<-filter(df,IL6_status == "H")
  df_IL6_hi$IL2_status <- ''
  df_IL6_hi[df_IL6_hi$IL2 > median(df_IL6_hi$IL2),]$IL2_status <- "H"
  df_IL6_hi[df_IL6_hi$IL2 <= median(df_IL6_hi$IL2),]$IL2_status <- "L"
  
  df_IL6_lo<-filter(df,IL6_status == "L")
  df_IL6_lo$IL2_status <- ''
  df_IL6_lo[df_IL6_lo$IL2 > median(df_IL6_lo$IL2),]$IL2_status <- "H"
  df_IL6_lo[df_IL6_lo$IL2 <= median(df_IL6_lo$IL2),]$IL2_status <- "L"
  
 TCGAanalyze_survival(df_IL2_hi,
                       clusterCol="IL6_status",
                       risk.table = FALSE,
                       conf.int = FALSE,
                       color = c("Dark2"),
                       filename = paste("TCGA",i,"IL2_hi.pdf",sep = "_"))  

 TCGAanalyze_survival(df_IL2_lo,
                      clusterCol="IL6_status",
                      risk.table = FALSE,
                      conf.int = FALSE,
                      color = c("Dark2"),
                      filename = paste("TCGA",i,"IL2_lo.pdf",sep = "_"))
 
 
 TCGAanalyze_survival(df_IL6_hi,
                      clusterCol="IL2_status",
                      risk.table = FALSE,
                      conf.int = FALSE,
                      color = c("Dark2"),
                      filename = paste("TCGA",i,"IL6_hi.pdf",sep = "_"))
 
 TCGAanalyze_survival(df_IL6_lo,
                      clusterCol="IL2_status",
                      risk.table = FALSE,
                      conf.int = FALSE,
                      color = c("Dark2"),
                      filename = paste("TCGA",i,"IL6_lo.pdf",sep = "_"))
 
 
 
 }



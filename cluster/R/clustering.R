

###
rendersilhouette <- function(data =  mtcars, k_start=2,k_end=5,
                                   pca_factors=1:ncol(data))
  {
  #pca on what variables
  data.pca <- prcomp(data[,pca_factors], center = TRUE,scale. = TRUE)
  
  multi_k = data.frame(k=k_start:k_end,mean_sil=0)
  
  
  for(k in k_start:k_end){
    km <- kmeans(data.pca$x[,1:10], centers = k, nstart=25)
    ss <- silhouette(km$cluster, dist(data.pca$x))
    multi_k[k-1,2] = mean(ss[, 3])
  }
  plot(multi_k[,1], type='b', multi_k[,2],main="Silhouette score \n Mesotheliama gene expression dataset", xlab='Number of clusters (k)', ylab='Average Silhouette Scores', frame=FALSE)
}

###
renderpca = function(data = meso[2],pcX=1,pcY=2,k = 2,
                    pca_factors=1:ncol(data))
  {
  #pca on what variables
  data.pca <- prcomp(data[,pca_factors], center = TRUE,scale. = TRUE)

  km <- kmeans(data.pca$x[,1:10], centers = k, nstart=25)
  
  colors = rainbow(length(unique(km$cluster)))
  
  colors_val <- lapply(data.frame(km$cluster), function(x) as.character(factor(km$cluster, labels=colors)))
  colors_val = as.character(colors_val$km.cluster)
  
  plot(data.pca$x[,pcX],data.pca$x[,pcY],xlab = paste0("PC",pcX),ylab = paste0("PC",pcY),col=colors_val,pch=8,main = "PCA & k grouping \n mesotheliama gene expression dataset")
  #biplot(data,choices=c(pcX,pcY))
  }



read_mesothelioma = function(path = "data/"){
###meso data:https://cran.r-project.org/web/packages/dnapath/vignettes/package_data.html#meso-data
  
  
  if(file.exists(paste0(path,"rnaseq_most_expr_meso.csv")) == F){
    library(readr)

    file_clinical <- paste0("http://linkedomics.org/data_download/TCGA-MESO/Human__",
                        "TCGA_MESO__MS__Clinical__Clinical__01_28_2016__BI__",
                        "Clinical__Firehose.tsi")
    file_rnaseq <- paste0("http://linkedomics.org/data_download/TCGA-MESO/Human__",
                      "TCGA_MESO__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene",
                      "__Firehose_RSEM_log2.cct.gz")

    clinical <- readr::read_tsv(file_clinical)
    rnaseq <- readr::read_tsv(file_rnaseq)
    write.csv(clinical,paste0(path,"clinical_meso.csv"),row.names=F)
  
    ###most expressed, 10% most variables genes (smaller dataset, keeping a lot of the signal)
    rnaseq_most_expressed = as.data.frame(rnaseq)
    rnaseq_most_expressed = rnaseq_most_expressed[rowMeans(rnaseq_most_expressed[,-1])>5.2,]
  
    sd = rep(0,nrow(rnaseq_most_expressed))
    for(i in 1:nrow(rnaseq_most_expressed))
      {
      sd[i] = sd(rnaseq_most_expressed[i,-1])
    }
    rnaseq_most_expressed_var = rnaseq_most_expressed[sd>quantile(sd,seq(0,1,by=0.1))[10],]
    write.csv(rnaseq_most_expressed_var,paste0(path,"rnaseq_most_expr_meso.csv"),row.names=F)
  }
  
  if(file.exists(paste0(path,"rnaseq_most_expr_meso.csv")) == T)
  
  rnaseq_most_expressed = read.csv(paste0(path,"rnaseq_most_expr_meso.csv"),row.names = 1)
  clinical = read.csv(paste0(path,"clinical_meso.csv"),row.names = 1)
  clinical = t(clinical)
  clinical = clinical[,colnames(clinical) %in% c("years_to_birth","pathologic_stage","histological_type","gender","overall_survival")]
  clinical = as.data.frame(clinical)
  clinical$years_to_birth = as.numeric(clinical$years_to_birth)
  clinical$overall_survival = as.numeric(clinical$overall_survival)
  
  list(clinical,t(rnaseq_most_expressed))
}


render_summary_data = function(clinical = clinical,variable = colnames(clinical)) {
  if(variable == c("overall_survival"))
  {
    hist(clinical[,variable],breaks =10,xlab = 'days of survival',ylab = "Number of patients",main = "Survival rates since diagnosis")
  }
  
  if(variable == c("years_to_birth"))
  {
    hist(clinical[,variable],breaks =10,xlab = 'Age (years)',ylab = "Number of patients",main = "Age at diagnosis")
  }
  
  if(variable == c("pathologic_stage"))
  {
    summary = rle(sort(clinical[,variable]))
    summary = data.frame(Pathologic_stage = summary$values,Nb_of_patients = summary$lengths)
    barplot(Nb_of_patients ~ Pathologic_stage, data = summary,main ="pathologic_stage")
  }
  
  if(variable == c("histological_type"))
  {
    summary = rle(sort(clinical[,variable]))
    summary = data.frame(histological_type = summary$values,Nb_of_patients = summary$lengths)
    barplot(Nb_of_patients ~ histological_type, data = summary,main ="histological_type")
    
  }
  
  if(variable == c("gender"))
  {
    summary = rle(sort(clinical[,variable]))
    summary = data.frame(gender = summary$values,Nb_of_patients = summary$lengths)
    barplot(Nb_of_patients ~ gender, data = summary,main ="gender")
    
  }
  
  
  
}


library(cluster)
library(DT)
library(wesanderson)



wes_colors = c(wes_palettes$GrandBudapest1,wes_palettes$GrandBudapest2,wes_palettes$Zissou1,wes_palettes$Rushmore)
###
rendersilhouette <- function(data.pca =  data.pca, k_start=2,k_end=5,
                                   pca_factors=1:ncol(data), colors = wes_colors)
  {
  #k results
  multi_k = data.frame(k=k_start:k_end,mean_sil=0)
  
  
  for(k in k_start:k_end){
    km <- kmeans(data.pca$x[,1:10], centers = k, nstart=25)
    ss <- silhouette(km$cluster, dist(data.pca$x))
    multi_k[k-1,2] = mean(ss[, 3])
  }
  plot(multi_k[,1], type='b', col=colors, pch = 19, lwd = 3,multi_k[,2],main="Silhouette score \n gene expression dataset", xlab='Number of clusters (k)', ylab='Average Silhouette Scores', frame=FALSE)
}

###
renderpca = function(data.pca=data.pca,pcX=1,pcY=2,k = 2,colors = wes_colors)
  {
  #kmeans
  km <- kmeans(data.pca$x[,1:10], centers = k, nstart=25)
  colors= colors[1:length(unique(km$cluster))]
  
  colors_val <- lapply(data.frame(km$cluster), function(x) as.character(factor(km$cluster, labels=colors)))
  colors_val = as.character(colors_val$km.cluster)
  
  plot(data.pca$x[,pcX],data.pca$x[,pcY],xlab = paste0("PC",pcX),ylab = paste0("PC",pcY),col=colors_val,pch=19,lwd=3,main = "PCA & k grouping \n gene expression dataset")
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
  
    #final cleanup
    rnaseq_most_expressed = read.csv(paste0(path,"rnaseq_most_expr_meso.csv"),row.names = 1)
    clinical = read.csv(paste0(path,"clinical_meso.csv"),row.names = 1)
    clinical = t(clinical)
    clinical = clinical[,colnames(clinical) %in% c("years_to_birth","pathologic_stage","histological_type","gender","overall_survival")]
    clinical = as.data.frame(clinical)
    clinical$years_to_birth = as.numeric(clinical$years_to_birth)
    clinical$overall_survival = as.numeric(clinical$overall_survival)
    clinical = clinical[match(colnames(rnaseq_most_expressed),rownames(clinical)),]
    colnames(clinical) = c("age","stage","histology","gender","survival")
    write.csv(clinical,paste0(path,"clinical_meso.csv"),row.names=T)
    
    }
  
  if((file.exists(paste0(path,"rnaseq_most_expr_meso.csv")) == T)) {
    rnaseq_most_expressed = read.csv(paste0(path,"rnaseq_most_expr_meso.csv"),row.names = 1)
    clinical = read.csv(paste0(path,"clinical_meso.csv"),row.names = 1)
    
  
  temp = list(clinical,t(rnaseq_most_expressed))}

  
}

read_data_gexp = function(file=NULL){
  
  if(!is.null(file)) {
    ext <- tools::file_ext(file)
    
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    temp = t(read.csv(file, header = T,row.names = 1))
  }
  temp
}



read_data_clinical = function(file=NULL){
  
  if(!is.null(file)) {
    ext <- tools::file_ext(file)
    
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    temp = read.csv(file, header = T,row.names = 1)
  }
  temp
}
  
  


render_summary_data = function(clinical = clinical,variable = colnames(clinical),colors = wes_colors) {
  if(variable == c("survival"))
  {
    hist(clinical[,variable],breaks =10,col = colors,xlab = 'days of survival',ylab = "Number of patients",main = "Survival rates since diagnosis")
  }
  
  if(variable == c("age"))
  {
    hist(clinical[,variable],breaks =10,col = colors,xlab = 'Age (years)',ylab = "Number of patients",main = "Age at diagnosis")
  }
  
  if(variable == c("stage"))
  {
    summary = rle(sort(clinical[,variable]))
    summary = data.frame(stage = summary$values,Nb_of_patients = summary$lengths)
    barplot(Nb_of_patients ~ stage, data = summary,main ="stage",col = colors,ylab = "Number of patients")
  }
  
  if(variable == c("histology"))
  {
    summary = rle(sort(clinical[,variable]))
    summary = data.frame(histology = summary$values,Nb_of_patients = summary$lengths)
    barplot(Nb_of_patients ~ histology, data = summary,main ="histology",col = colors,ylab = "Number of patients")
    
  }
  
  if(variable == c("gender"))
  {
    summary = rle(sort(clinical[,variable]))
    summary = data.frame(gender = summary$values,Nb_of_patients = summary$lengths)
    barplot(Nb_of_patients ~ gender, data = summary,main ="gender",col = colors,ylab = "Number of patients")
    
  }
  
  if(variable == c("cluster"))
  {
    summary = rle(sort(clinical[,variable]))
    summary = data.frame(cluster = summary$values,Nb_of_patients = summary$lengths)
    barplot(Nb_of_patients ~ cluster, data = summary,main ="cluster",col = colors,ylab = "Number of patients")
    
  }
}


render_detailed_silhouette = function(input_pca_data = data.pca,k = 2,colors = wes_colors){
  km <- kmeans(input_pca_data$x[,1:10], centers = k, nstart=25)
  ss <- silhouette(km$cluster, dist(input_pca_data$x))
  plot(ss,col = colors[1:k],nmax = 10,main = paste0("Detailed Silhouette plot for k = ",k))
}




toto = function(data.pca=data.pca,k = 2,colors = wes_colors){
  
  plot(1:10)
}
  
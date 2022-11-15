
library(cluster)
library(DT)
library(wesanderson)



wes_colors = c(wes_palettes$GrandBudapest1,wes_palettes$GrandBudapest2,wes_palettes$Zissou1,wes_palettes$Rushmore)
###


clustering_statistics <- function(input_data =  data[[2]], k_start=2,k_end=5,type= "kmeans")
{
  #k results
  input_pca_data.pca = prcomp(input_data[,1:ncol(input_data)], center = TRUE,scale. = TRUE)
  multi_k = data.frame(k=k_start:k_end,mean_sil=0,twss=0)
  
  if(type == "kmeans"){
    
    for(k in k_start:k_end){
      km <- kmeans(input_pca_data.pca$x[,1:10], centers = k, nstart=25)
      ss <- silhouette(km$cluster, dist(input_pca_data.pca$x))
      multi_k[k-1,2] = mean(ss[, 3])
      multi_k[k-1,3] = km$tot.withinss
    }
  }
  
  if(type == "dendrogram"){
    for(k in k_start:k_end){
      cl <- hclust(dist(input_data))
      km <- list(cluster = cutree(cl, k=k))
      ss <- silhouette(cutree(cl, k=k) ,dist(input_data), title=title(main = 'Good'))
      multi_k[k-1,2] = mean(ss[, 3])
      
      #this is a crutch for now...
      multi_k[k-1,3] = kmeans(input_pca_data.pca$x[,1:10], centers = k, nstart=25)$tot.withinss
    }
  }
  multi_k
}
  
  
#####  
renderstwss <- function(multi_k =  clustering_statistics(data[[2]]),
                        colors = wes_colors)
{
  plot(multi_k[,1],multi_k[,3], type='b', col=colors, pch = 19, lwd = 3,main="Total within sum of square \n gene expression dataset", xlab='Number of clusters (k)', ylab='Total Within Sum of Square', frame=FALSE)
}


rendersilhouette <- function(multi_k =  clustering_statistics(data[[2]]),
                             colors = wes_colors)
  {
  plot(multi_k[,1],multi_k[,2], type='b', col=colors, pch = 19, lwd = 3,main="Silhouette score \n gene expression dataset", xlab='Number of clusters (k)', ylab='Average Silhouette Scores', frame=FALSE)
}



render_detailed_silhouette = function(input_data = data[[2]],k = 2,colors = wes_colors,type=c("kmeans","dendrogram")){
  if(type == "kmeans"){
    input_pca_data.pca = prcomp(input_data[,1:ncol(input_data)], center = TRUE,scale. = TRUE)
    km <- kmeans(input_pca_data.pca$x[,1:10], centers = k, nstart=25)
    ss <- silhouette(km$cluster, dist(input_pca_data.pca$x))
    rownames(ss) = names(km$cluster)
  }
  
  if(type == "dendrogram"){
    cl <- hclust(dist(input_data))
    ss <- silhouette(cutree(cl, k=k) ,dist(input_data), title=title(main = 'Good'))
    rownames(ss) = cl$labels
  }
  plot(ss,col = colors[1:k],max.strlen=20,nmax.lab = 200,cex.names = 0.2,main = paste0("Detailed Silhouette plot for k = ",k)) 
}


###
renderpca = function(input_data=data[[2]],pcX=1,pcY=2,k = 2,colors = wes_colors,type="kmeans")
  {
  #kmeans
  input_pca_data.pca = prcomp(input_data[,1:ncol(input_data)], center = TRUE,scale. = TRUE)
  if(type == "kmeans"){
    km <- kmeans(input_pca_data.pca$x[,1:10], centers = k, nstart=25)
  }
  if(type == "dendrogram"){
    cl <- hclust(dist(input_data))
    km <- list(cluster = cutree(cl, k=k))
  }
  
  
  colors= colors[1:length(unique(km$cluster))]
  
  colors_val <- lapply(data.frame(km$cluster), function(x) as.character(factor(km$cluster, labels=colors)))
  colors_val = as.character(colors_val$km.cluster)
  eigen_vector = input_pca_data.pca$sdev^2 
  pve= round(eigen_vector/sum(eigen_vector)*100,2)
  plot(input_pca_data.pca$x[,pcX],input_pca_data.pca$x[,pcY],xlab = paste0("PC",pcX," (",pve[pcX]," %)"),ylab =paste0("PC",pcY," (",pve[pcY]," %)"),col=colors_val,pch=19,lwd=3,main = paste0("Principal Component Analysis \n (grouping according to ",type,")"))
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
    rnaseq_most_expressed_var = rnaseq_most_expressed[sd>quantile(sd,seq(0,1,by=0.1))[5],]
    write.csv(rnaseq_most_expressed_var,paste0(path,"rnaseq_most_expr_meso50.csv"),row.names=F)
  
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
  
  if(variable == c("kmeans"))
  {
    summary = rle(sort(clinical[,variable]))
    summary = data.frame(kmeans = summary$values,Nb_of_patients = summary$lengths)
    barplot(Nb_of_patients ~ kmeans, data = summary,main ="Clusters kmeans",col = colors,ylab = "Number of patients")
    
  }
  
  if(variable == c("hierarchical"))
  {
    summary = rle(sort(clinical[,variable]))
    summary = data.frame(hierarchical = summary$values,Nb_of_patients = summary$lengths)
    barplot(Nb_of_patients ~ hierarchical, data = summary,main ="Clusters hierarchical",col = colors,ylab = "Number of patients")
    
  }
}


hc = function(gexp = data[[2]],max_genes=200){
  gexp_max = t(gexp)[1:max_genes,]
  stats::heatmap(gexp_max,main=list("Heatmap (genes X samples)",cex = 1.1))
  }


plots = function(multi_k = multi_k, data = data,clustering_metrics = '',k = 2) {
  switch(clustering_metrics,
       silhouette = rendersilhouette(multi_k = multi_k),
       detailed = render_detailed_silhouette(input_data = data[[2]],k = k,type = 'kmeans'),
       twss = renderstwss(multi_k = multi_k),
       silhouette_hc = rendersilhouette(multi_k = multi_k),
       detailed_hc = render_detailed_silhouette(input_data = data[[2]],k = k,type = 'dendrogram'),
       twss_hc = renderstwss(multi_k = multi_k)
  )

}


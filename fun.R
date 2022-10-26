
#fun1
renderphylo = function(drop="",type = "phylogram"){
  #Slightly better than default plot
  if(drop!="DON'T DROP") bird.orders_dropped= drop.tip(bird.orders,drop)
  if(drop=="DON'T DROP") bird.orders_dropped= bird.orders
  if(drop=="DROP RANDOMLY 10 TIPS") bird.orders_dropped= drop.tip(bird.orders,bird.orders$tip.label[sample(1:23,10)])
  phylo2 = plot.phylo(bird.orders_dropped,type=type, edge.width = 1, label.offset = 1,cex = 0.4, main = "Figure 1 (dropping tips/type)")
  
}


#fun2
rendercoloured = function(color1="blue", color2="red",pies = "YES",tips_colors= "YES") {
    trait1 = c(sample(c(0:25),6),c(sample(25:75,6)),sample(c(75:100),11))
    names(trait1) = bird.orders$tip.label
    
    #Use contMap from phytools
    tree.trait <- contMap(bird.orders, trait1, plot = F)
    tree.trait<-setMap(tree.trait, colors=c(color1,color2))
    
    #Carefull, we now plot using the plotting tool from phytools
    plot.contMap(tree.trait, type = "phylogram", legend = 0.4*max(nodeHeights(bird.orders)), fsize = c(0.9, 0.7), lwd=3, edge.width = 2, offset = 1,mar = c(2,4,3,2))
    
    #add pies
    if(pies=="YES"){
      pies = matrix(0,nrow =bird.orders$Nnode,ncol=2)
      for(i in 1:bird.orders$Nnode) {x = sample(seq(0,1,by = 0.001),1); pies[i,] = c(x,1-x)}
      nodelabels(pie = pies, piecol = c(color1,color2), cex = 0.6)
    }
    
    if(tips_colors=="YES"){  
      tip.labels.col = c(rep(color1,11),rep(color2,12))
      tiplabels(pch=15, cex=1.2, col= tip.labels.col,offset = 0.7)
    }
    
    #Add a title manually to plot.contMap
    title("Figure 2 (adding branch colors according to a trait)",xpd = T)
  
  
}


#fun3
rendermap = function() {

renderPlot({
  
  #Matrix of random USA GPS coordinates 
  mat = matrix(c(sample(seq(from = 32, to = 45, by = 0.001),23),sample(seq(from = -120, to = -84, by = 0.001),23)),nrow = 23,ncol =2)
  rownames(mat) = bird.orders$tip.label
  
  ###World mapping function (Here we are using a USA map, but this can be changes with the database option)
  phylo.to.map(bird.orders,coords = mat, database = "usa")
  
  #Add a title manually to plot.contMap
  title("Figure 4 (geolocalised samples map)",line = -3)
},width=1400, height = 1200)

}


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
  plot(multi_k[,1], type='b', multi_k[,2],main="Clustering (silhouette score) \n Mesotheliama gene expression dataset", xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)
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



read_mesothelioma = function(path = "~/Documents/git_repos/shiny_examples/data/"){
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
  list(clinical,t(rnaseq_most_expressed))
}

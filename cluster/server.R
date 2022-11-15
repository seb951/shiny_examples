#
library(shiny)
library(cluster)
library(DT)
library(wesanderson)
source("R/clustering.R")

shinyServer(function(input, output,session) {
  #get the data ready
    observe({
      if(is.null(input$GEXP_file)){
        data = read_mesothelioma()}
      
      if(!is.null(input$GEXP_file) & is.null(input$CLINICAL_file)){
        data_gexp = read_data_gexp(file=input$GEXP_file$datapath)
        data_clinical = data.frame(age = rep(0,nrow(data_gexp)), stage = 0,histology = 0, gender = 0,survival = 0)
        data = list(data_clinical,data_gexp)}
      
      if(!is.null(input$GEXP_file) & !is.null(input$CLINICAL_file)) {
        data_gexp = read_data_gexp(file=input$GEXP_file$datapath)
        data_clinical = read_data_clinical(file=input$CLINICAL_file$datapath)
        data = list(data_clinical,data_gexp)                     
        }
 
        
    data.pca <- prcomp(data[[2]][,1:ncol(data[[2]])], center = TRUE,scale. = TRUE)
    data.dist = dist(data[[2]])
    cl <- hclust(data.dist)
    
    
    #make plots
    output$pca = renderPlot({
        renderpca(data[[2]],
                 pcX=as.numeric(input$pc[1]),
                 pcY=as.numeric(input$pc[2]),
                 k = input$k_selected,
                 type= "kmeans")
    })
    
    output$pca_hc = renderPlot({
      renderpca(data[[2]],
                pcX=1,
                pcY=2,
                k = input$k_selected_hc,
                type = "dendrogram")
    })
    
    multi_k = clustering_statistics(data[[2]],
                                    k_start=2,
                                    k_end=10,
                                    type = "kmeans")
    
    multi_k_hc = clustering_statistics(data[[2]],
                                       k_start=2,
                                       k_end=10,
                                       type = "dendrogram")
    
    if(input$clustering_metrics == 'silhouette') {
      output$clustering_plots = renderPlot({rendersilhouette(multi_k = multi_k)})
    }
    
    if(input$clustering_metrics == 'twss') {
      output$clustering_plots = renderPlot({renderstwss(multi_k = multi_k_hc)})
    }
    
    if(input$clustering_metrics == 'detailed') {
      output$clustering_plots = renderPlot({render_detailed_silhouette(input_data = data[[2]],k = input$k_selected,type = "kmeans")})
    }
    
    ####hierarchical plots
    if(input$hc_metrics == 'silhouette_hc') {
      output$hc_plots = renderPlot({rendersilhouette(multi_k = multi_k_hc)})
    }
    
    if(input$hc_metrics == 'twss_hc') {
      output$hc_plots = renderPlot({renderstwss(multi_k = multi_k_hc)})
    }
    
    if(input$hc_metrics == 'detailed_hc') {
      output$hc_plots = renderPlot({render_detailed_silhouette(input_data = data[[2]],k = input$k_selected_hc,type = "dendrogram")})
    }

    output$heatmap = renderPlot({
      hc(gexp = data[[2]])})
    
    observe({
        if(length(input$pc) > 2){
            selected = tail(input$pc,2)
            updateCheckboxGroupInput(session, "pc", selected= selected)
            print(selected)
        }
        if(length(input$pc) == 1){
            selected = c(1,ifelse(input$pc==1,2,1))
            updateCheckboxGroupInput(session, "pc", selected = selected)
        }
    })
    
    output$downloadGexp <- downloadHandler(
      filename = function() {
        paste('Gexpr_data_', Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        write.csv(data[[2]], con)
      }
    ) 
    
    output$downloadHeatmap <- downloadHandler(
      filename = function() { paste0("heatmap_",Sys.Date(),'.png') },
      content = function(file) {
        png(file,res= 200,width = 800, height = 800)
        hc(gexp = data[[2]]) 
        dev.off()
      }
    ) 
    
    
    
    ####CLINICAL
    observe({
      km = kmeans(data.pca$x[,1:10], centers = input$k_selected, nstart=25)
      sil_cl <- silhouette(cutree(cl, k=input$k_selected_hc) ,data.dist, title=title(main = 'Good'))
      clinical2 = data.frame(kmeans=km$cluster,hierarchical=sil_cl[,2],data[[1]])
      output$mesotable = DT::renderDataTable({clinical2})
      
      output$summary_data = renderPlot({
        render_summary_data(clinical = clinical2,variable = input$variable) })
        
      
        output$downloadData <- downloadHandler(
          filename = function() {
            paste('clinical_cluster_data_', Sys.Date(), '.csv', sep='')
          },
          content = function(con) {
            write.csv(clinical2, con)
          }
        )    
      
    })
    
        
})
    
    ###TEXT
    output$general <- renderUI({
      para7 <- "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Here, I use a publicly available <a href='https://cran.r-project.org/web/packages/dnapath/vignettes/package_data.html#meso-data'> dataset</a> as an example on how you can cluster gene expression data and report results interactively.
      This dataset contains gene expression and clinical data obtained from 87 cancer (mesothelioma) patients. <a href='https://en.wikipedia.org/wiki/Mesothelioma'>Mesothelioma </a> is a type of cancer that develops from the thin layer of tissue that covers many of the internal organs (the mesothelium) and mainly affects lungs.
      <br/><br/>
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Alternatively, you can upload your own gene expression dataset (<i>.csv</i> file of normalized counts) and clinical data (<i>.csv</i> file), using the input buttons to the left.
      <br/><br/>"
      HTML(para7)
    })
    
    output$kmeans <- renderUI({
      para7 <- "Below are the results of clustering based on kmeans, silhouette score & TWSS to choose an optimal k value and pca for visualisation. I use a (VST normalised) gene expression dataset, where I keep only the 10% most variables genes above a specific expression threshold (4.2).
      <br/><br/>" 
      HTML(para7)
    })
    
    output$dendrogram <- renderUI({
      para7 <- "Below are the results of clustering based on cuting a sample dendrogram in a fixed number of groups k, silhouette score & TWSS to choose an optimal k value and pca for visualisation. I use a (VST normalised) gene expression dataset, where I keep only the 10% most variables genes above a specific expression threshold (4.2). 
      Clustering based .
      <br/><br/>"
      HTML(para7)
    })
    
    output$clinical <- renderUI({
        para7 <- "Below are the clinical data from the mesothelioma example dataset, both in summarized plots and in a downloadable datatable (only a relevant subset of clinical data is used for simplicity). 
        <br/><br/> Note that these ouputs are reactive based on the K decisions in the two previous tabs.
        <br/><br/>"
        HTML(para7)
        })
    
    
    
    
    
})




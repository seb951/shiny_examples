#
library(shiny)
library(cluster)
library(DT)
library(wesanderson)
source("R/clustering.R")

shinyServer(function(input, output,session) {
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
    
    
    output$sil = renderPlot({
        rendersilhouette(data[[2]],
                         k_start=2,
                         k_end=10,
                         type = "kmeans")
    })
    
    output$sil_hc = renderPlot({
      rendersilhouette(data[[2]],
                       k_start=2,
                       k_end=10,
                       type = "dendrogram")
    })
    
    output$detailedsil = renderPlot({
      render_detailed_silhouette(input_data = data[[2]],k = input$k_selected,type = "kmeans")
    })
    
    output$heatmap = renderPlot({
      hc(gexp = data[[2]])})
    
    output$detailedsil_hc = renderPlot({
      render_detailed_silhouette(input_data = data[[2]],k = input$k_selected_hc,type = "dendrogram")})
    

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
    
    ####CLINICAL
    observe({
      km = kmeans(data.pca$x[,1:10], centers = input$k_selected, nstart=25)
      sil_cl <- silhouette(cutree(cl, k=input$k_selected_hc) ,data.dist, title=title(main = 'Good'))
      clinical2 = data.frame(kmeans=km$cluster,dendrogram=sil_cl[,2],data[[1]])
      output$mesotable = DT::renderDataTable({clinical2})
      
      output$summary_data = renderPlot({
        render_summary_data(clinical = clinical2,variable = input$variable) 
        
      
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
    
        
})
    
    ###TEXT
    output$general <- renderUI({
      para7 <- "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Here, I use a publicly available <a href='https://cran.r-project.org/web/packages/dnapath/vignettes/package_data.html#meso-data'> dataset</a> as an example on how you can cluster gene expression data and report results interactively.
      This dataset contains gene expression and clinical data obtained from 87 cancer (mesothelioma) patients. <a href='https://en.wikipedia.org/wiki/Mesothelioma'>Mesothelioma </a> is a type of cancer that develops from the thin layer of tissue that covers many of the internal organs (the mesothelium) and mainly affects lungs.
      <br/><br/>
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Alternatively, you can upload your own gene expression dataset (<b>genes X patients</b>, <i>.csv</i> file of normalized counts) and clinical data (<b>patients X attributes</b>, <i>.csv</i> file) in the Kmeans tab, using the input buttons to the left.
      <br/><br/>"
      HTML(para7)
    })
    
    output$kmeans <- renderUI({
      para7 <- "Below are the results of clustering based on kmeans,silhouette score to choose an optimal k value and pca for visualisation. I use a (VST normalised) gene expression dataset, where I keep only the 10% most variables genes above a specific expression threshold (4.2).
      <br/><br/>" 
      HTML(para7)
    })
    
    output$dendrogram <- renderUI({
      para7 <- "Clustering based on cuting the sample dendrogram in a fixed number of groups k. Note that heatmap contains only 200 genes for cleaner visualisation.
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

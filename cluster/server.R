#
library(shiny)
library(cluster)
library(DT)
library(wesanderson)
source("R/clustering.R")

shinyServer(function(input, output,session) {
    observe({
      if(is.null(input$GEXP_file)){
        data = read_mesothelioma()
        print("gexpNULL,clinicalEITHER")}
      if(!is.null(input$GEXP_file) & is.null(input$CLINICAL_file)){
        data_gexp = read_data_gexp(file=input$GEXP_file$datapath)
        data_clinical = data.frame(age = rep(0,nrow(data_gexp)), stage = 0,histology = 0, gender = 0,survival = 0)
        data = list(data_clinical,data_gexp)      
        print("gexpNOTNULL,clinicalNULL")}
      if(!is.null(input$GEXP_file) & !is.null(input$CLINICAL_file)) {
        print(input$GEXP_file)
        print(input$CLINICAL_file)
        data_gexp = read_data_gexp(file=input$GEXP_file$datapath)
        data_clinical = read_data_clinical(file=input$CLINICAL_file$datapath)
        data = list(data_clinical,data_gexp)                     
        print("NOTNULL")}
 
        
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
    output$credential <- renderUI({
        para4 <- "&nbsp;&nbsp;&nbsp;sebastien.renaut@gmail.com, 2022"
        HTML(para4,sep = '<br/><br/>')
    })
    
    
    output$kmeans <- renderUI({
      para7 <- "I use a (VST normalised) mesothelioma gene expression dataset from <a href='https://cran.r-project.org/web/packages/dnapath/vignettes/package_data.html#meso-data'>here</a>.
        I filtered out low expressed genes and kept only 10% most variable genes to speed things up. 
        But you can upload your own dataset using the input buttons to the left. Below are the results of clustering based on kmeans, pca and silhouette score."
      para5 <- "<br/><br/>" 
      HTML(paste(para7,para5,sep = '<br/><br/>'))
    })
    
    output$dendrogram <- renderUI({
      para7 <- "Here is a simple example of clustering based on hierchical clustering and cutting a dendrogram with a defined number of groups (k)."
      para5 <- "<br/><br/>" 
      HTML(paste(para7,para5,sep = '<br/><br/>'))
    })
    
    output$clinical <- renderUI({
        para7 <- "Below are the clinical data, both in summarized graphs and in a datatable. You can download the datatable. Note that these are reactive based on the K decisions in the previous tabs"
        para5 <- "<br/><br/>" 
        HTML(paste(para7,para5,sep = '<br/><br/>'))
        })
    
    
    
    
    
})

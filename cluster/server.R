#
library(shiny)
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
    output$pca = renderPlot({
        renderpca(data.pca,
                 pcX=as.numeric(input$pc[1]),
                 pcY=as.numeric(input$pc[2]),
                 k = input$k_selected)
    })
    
    
    output$sil = renderPlot({
        rendersilhouette(data.pca,
                         k_start=2,
                         k_end=10)
    })
    output$detailedsil = renderPlot({
      render_detailed_silhouette(input_pca_data = data.pca,k = input$k_selected)
    })
    
    

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
      clinical2 = data.frame(cluster=km$cluster,data[[1]])
      output$mesotable = DT::renderDataTable({clinical2})
      
      output$summary_data = renderPlot({
        render_summary_data(clinical = clinical2,variable = input$variable) 
        
        
      })
    })
    
        
})
    
    ###TEXT
    output$general <- renderUI({
        if(input$dataset_name != "mesothelioma") 
        {para4 <- paste0("Here is a simple example of clustering based on kmeans, pca and silhouette score. I use a user-inputed dataset called: ",input$dataset_name)}
        else {para4 <- "Here is a simple example of clustering based on kmeans, pca and silhouette score. I use a mesothelioma clinical and gene expression dataset from <a href='https://cran.r-project.org/web/packages/dnapath/vignettes/package_data.html#meso-data'>here</a>."}
        
        
        
        para5 <- "<br/><br/>" 
        
        HTML(paste(para4,para5, sep = '<br/><br/>'))
        
    })
    
    output$gene_expression <- renderUI({
        para5 <- "<br/><br/>" 
        if(input$dataset_name == "mesothelioma") 
        {para6 = "Below, I use VST normalised gene expression data. I filtered out low expressed genes and kept only 10% most variable genes to speed things up. Off course, it's easy to change these filter or make them reactive."
        HTML(paste(para6,para5, sep = '<br/><br/>'))}
        else {
            HTML(paste("", sep = '<br/><br/>'))
        }
        
    })
    
    
    
    
    output$clinical <- renderUI({
        if(input$dataset_name == "mesothelioma") 
        {para7 <- "Below I add a data table and summary plots. Note that these are reactive based on the K decision above"
        para5 <- "<br/><br/>" 
        HTML(paste(para7,para5sep = '<br/><br/>'))} else {
            HTML(paste("",sep = '<br/><br/>'))
        }
            
        
    })
    
    
    
    
    
})

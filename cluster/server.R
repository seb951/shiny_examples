#
library(shiny)
source("R/clustering.R")

shinyServer(function(input, output,session) {
    
    meso = read_mesothelioma()
    colnames(meso[[1]]) = c("age","stage","histology","gender","survival")
    data.pca <- prcomp(meso[[2]][,1:ncol(meso[[2]])], center = TRUE,scale. = TRUE)
    
    
    output$pca = renderPlot({
        renderpca(data.pca,
                 pcX=as.numeric(input$pc[1]),
                 pcY=as.numeric(input$pc[2]),
                 k = input$k_selected)
    })
    output$sil = renderPlot({
        rendersilhouette(data.pca,
                         k_start=input$krange[1],
                         k_end=input$krange[2])
    })

    observe({
        if(length(input$pc) > 2){
            updateCheckboxGroupInput(session, "pc", selected= tail(input$pc,2))
        }
        if(length(input$pc) == 1){
            updateCheckboxGroupInput(session, "pc", selected= c(1,ifelse(input$pc==1,2,1)))
        }
    })
        
    observe({
        km = kmeans(data.pca$x[,1:10], centers = input$k_selected, nstart=25)
        clinical = data.frame(cluster=km$cluster,meso[[1]])
        output$mesotable = DT::renderDataTable({clinical})
        
        output$summary_data = renderPlot({
            render_summary_data(clinical = clinical,variable = input$variable) 
            
        
        })
    
       
        
    
    
    
    ###
    output$general <- renderUI({
        para4 <- "Here is a simple example of clustering based on kmeans, pca and silhouette score. I use a mesothelioma clinical and gene expression dataset from <a href='https://cran.r-project.org/web/packages/dnapath/vignettes/package_data.html#meso-data'>here</a>."
        para5 <- "<br/><br/>" 
        
        HTML(paste(para4,para5, sep = '<br/><br/>'))
        
    })
    
    output$gene_expression <- renderUI({
        para6 = "Below, I use VST normalised gene expression data. I filtered out low expressed genes and kept only 10% most variable genes to speed things up. Off course, it's easy to change these filter or make them reactive."
        para5 <- "<br/><br/>" 
        
        HTML(paste(para6,para5, sep = '<br/><br/>'))
        
    })
    
    output$clinical <- renderUI({
        para7 <- "Below I add a data table and summary plots. Note that these are reactive based on the K decision above"
        para5 <- "<br/><br/>" 
        HTML(paste(para7,para5sep = '<br/><br/>'))
        
    })
    
    
     
        
        
        
    })
    
    
    
})

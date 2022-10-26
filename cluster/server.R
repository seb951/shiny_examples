#
library(shiny)
source("~/Documents/git_repos/shiny_examples/R/clustering.R")

shinyServer(function(input, output,session) {
    
    meso = read_mesothelioma()
                 
    output$mesotable = DT::renderDataTable({
        meso[[1]]})
    
    
    output$pca = renderPlot({
        renderpca(data=meso[[2]],
                 pcX=as.numeric(input$pc[1]),
                 pcY=as.numeric(input$pc[2]),
                 k = input$k_selected)
                 #pca_factors=input$pca_factors)
    })
    output$sil = renderPlot({
        rendersilhouette(data = meso[[2]],
                         k_start=input$krange[1],
                         k_end=input$krange[2])
                         #pca_factors=input$pca_factors)
    })

    observe({
        if(length(input$pc) > 2){
            updateCheckboxGroupInput(session, "pc", selected= tail(input$pc,2))
        }
        if(length(input$pc) == 1){
            updateCheckboxGroupInput(session, "pc", selected= c(1,ifelse(input$pc==1,2,1)))
        }
    })
    
    
    ###
    output$text_cluster <- renderUI({
        para4 <- "Below, I have a little example of clustering based on kmeans, pca and silhouette score"
        para5 <- "I use a mesothelimia gene expression dataset from <a href='https://cran.r-project.org/web/packages/dnapath/vignettes/package_data.html#meso-data'>here</a>. This contains VST normalised gene expression data. I filtered out low expressed genes and kept only 10% most variable genes to speed things up."
        
        HTML(paste(para4,para5, sep = '<br/><br/>'))
        
    })
    
    output$summary_data = renderPlot({
        render_summary_data(clinical = meso[[1]],variable = input$variable) 
    })
    
    
    
})

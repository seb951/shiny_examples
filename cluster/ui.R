#
library(shiny)
library(cluster)
library(DT)
library(wesanderson)
library(shinythemes)

shinyUI(navbarPage(title = "Clustering example",theme = shinytheme("cerulean"),
                   header = list(
    htmlOutput("credential")),
    tabPanel("Kmeans",
             sidebarPanel(
               fileInput("GEXP_file","Gene Expression dataset (.csv with header)", accept=".csv", width = NULL, placeholder = ""),
               fileInput("CLINICAL_file","Clinical dataset (.csv with header)", accept=".csv", width = NULL, placeholder = ""),
             ),
             mainPanel(h4("The gene expression dataset (Kmeans)"),
                       htmlOutput("kmeans"),
                       plotOutput("sil")),
             
             sidebarPanel(
                 sliderInput("k_selected","Select a K based on Silhouette plots",min = 2, max = 8, value = 2),
                 checkboxGroupInput("pc", "Choose two PCs for visualisation",
                                    choiceNames = paste0("PC",1:8),
                                    selected = c(1,2),
                                    choiceValues =  1:8)
             ),
             mainPanel(plotOutput("detailedsil"),
                       plotOutput("pca"))
    ),
    
    tabPanel("Hierarchical clustering",
             sidebarPanel(strong("Note that heatmap contains only 200 genes for cleaner visualisation")
             ),
             mainPanel(
                 h4("The gene expression dataset (Hierchical Clustering)"),
               htmlOutput("dendrogram"),
               
               plotOutput("heatmap"),
               
             ),
             sidebarPanel(sliderInput("k_selected_hc","Select a K to cut (samples) dendrogram ",min = 2, max = 8, value = 2)
                          ),
             mainPanel(plotOutput("sil_hc"),
                       plotOutput("detailedsil_hc"),
                       plotOutput("pca_hc"))
    ),
    
    
    
    
    tabPanel("Summary",
             sidebarPanel(
               radioButtons("variable", "Choose a clinical variable to display:",
                            selected = "age",
                            choices= c("kmeans","dendrogram","age","stage","histology","gender","survival")),
               downloadButton('downloadData', 'Download Clinical Data')
             ),
             
                mainPanel(
                  h4("The Clinical dataset (summary stats)"),
                 htmlOutput("clinical"),
                   plotOutput("summary_data")),
             mainPanel(
                
                 DT::dataTableOutput("mesotable"))
    )
    
    
    
)
)
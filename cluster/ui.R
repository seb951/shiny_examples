#
library(shiny)
library(cluster)
library(DT)
library(wesanderson)
library(shinythemes)

shinyUI(navbarPage(title = strong("Data Report"),
                   theme = shinytheme("cerulean"),
                   tags$head(
                       tags$style(type = "text/css", ".container-fluid {padding-left:10px;
                padding-right:10px; margin-right:10px; margin-left:10px;}")),
                   header = list(strong("General description"),htmlOutput("general")),
    tabPanel("Kmeans",
             sidebarPanel(
               fileInput("GEXP_file","Gene Expression dataset (.csv with header)", accept=".csv", width = NULL, placeholder = ""),
               fileInput("CLINICAL_file","Clinical dataset (.csv with header)", accept=".csv", width = NULL, placeholder = ""),
             ),
             mainPanel(h4("Gene expression clustering (Kmeans)"),
                       htmlOutput("kmeans"),
                       plotOutput("sil")),
             
             sidebarPanel(
                 sliderInput("k_selected","Select a K based on Silhouette plots",min = 2, max = 8, value = 2)
             ),
             mainPanel(plotOutput("detailedsil")),
             
             sidebarPanel(
                 checkboxGroupInput("pc", "Choose two PCs for visualisation",
                                    choiceNames = paste0("PC",1:8),
                                    selected = c(1,2),
                                    choiceValues =  1:8)
             ),
             mainPanel(plotOutput("pca")),
             tags$footer("sebastien.renaut@gmail.com --- 2022")
    ),
    
    tabPanel("Hierarchical clustering",
             mainPanel(width = 12,
                 h4("Gene expression clustering (Hierchical)"),
               htmlOutput("dendrogram"),
               
               plotOutput("heatmap"),
               
             ),
             sidebarPanel(sliderInput("k_selected_hc","Select a K based on Silhouette plots to cut dendrogram ",min = 2, max = 8, value = 2)
                          ),
             mainPanel(plotOutput("sil_hc"),
                       plotOutput("detailedsil_hc"),
                       plotOutput("pca_hc")),
             tags$footer("sebastien.renaut@gmail.com --- 2022")
    ),
    
    
    
    
    tabPanel("Summary",
             sidebarPanel(
                 radioButtons("variable", "Choose a clinical variable to display:",
                              selected = "age",
                              choices= c("kmeans","dendrogram","age","stage","histology","gender","survival")),
                 
                 downloadButton('downloadData', 'Download Clinical Data')
             ),
             mainPanel(
                 h4("The Clinical dataset"),
                 htmlOutput("clinical"),
                 plotOutput("summary_data")),
             mainPanel(
                 DT::dataTableOutput("mesotable",width = "66%"),tags$table(tableFooter("sebastien.renaut@gmail.com --- 2022")))
             
             
                 
    )
    
    
    
)
)
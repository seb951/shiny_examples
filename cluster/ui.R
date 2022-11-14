#
library(shiny)
library(cluster)
library(DT)
library(wesanderson)
library(shinythemes)

shinyUI(navbarPage(title = strong("Data Report"),
                   theme = shinytheme("cerulean"),
                   tabPanel("Summary",
                            icon = icon("clipboard"),
                            mainPanel(width =12,strong("General description"),
                                     htmlOutput("general")),
                            sidebarPanel(
                                fileInput("GEXP_file","Gene Expression dataset", accept=".csv", width = NULL, placeholder = "Gexp.csv"),
                                fileInput("CLINICAL_file","Clinical dataset", accept=".csv", width = NULL, placeholder = "Clinical.csv")),
                            mainPanel(
                                plotOutput("heatmap")),
                            tags$footer("sebastien.renaut@gmail.com --- 2022")),
                  
                   navbarMenu("Clustering",
                              icon = icon("circle-nodes"),
    tabPanel("K-means",
             mainPanel(width =12,h4("Gene expression clustering (K-means)"),
                      htmlOutput("kmeans")),
             
             mainPanel("",width =12,
                           fluidRow(
                               splitLayout(cellWidths = c("50%", "50%"), plotOutput("sil"), plotOutput("twss"))
                           )), 
             sidebarPanel("",
                 sliderInput("k_selected","Select a K for kmeans clustering",min = 2, max = 8, value = 2)),
             mainPanel("",
                 plotOutput("detailedsil")
             ),
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
             mainPanel(width =12,
                 h4("Gene expression clustering (Hierarchical)"),
               htmlOutput("dendrogram")
             ),
             mainPanel("",width =12,
                       fluidRow(
                           splitLayout(cellWidths = c("50%", "50%"), plotOutput("sil_hc"), plotOutput("twss_hc"))
                       )), 
             sidebarPanel(sliderInput("k_selected_hc","Select a K for Hierarchical clustering",min = 2, max = 8, value = 2)
                          ),
             mainPanel(plotOutput("detailedsil_hc"),
                       plotOutput("pca_hc")),
             tags$footer("sebastien.renaut@gmail.com --- 2022")
    )
                   ),
    
    
    
    tabPanel("Clinical",
             icon = icon("microscope"),
             sidebarPanel(
                 radioButtons("variable", "Choose a clinical variable to display:",
                              selected = "age",
                              choices= c("kmeans","hierarchical","age","stage","histology","gender","survival"))
                 
                
             ),
             mainPanel(
                 h4("The Clinical dataset"),
                 htmlOutput("clinical"),
                 plotOutput("summary_data")),
             mainPanel(
                 DT::dataTableOutput("mesotable",width = "66%"),
                 tags$table(tableFooter("sebastien.renaut@gmail.com --- 2022")))
    ),
    navbarMenu("Downloads",icon = icon("floppy-disk"),
    tabPanel("Clinical/Clustering",
             sidebarPanel(
                 downloadButton('downloadData', 'Download Clinical/Clustering'))),
    tabPanel("Gene Expression",
             sidebarPanel(
                 downloadButton('downloadGexp', 'Download Gene Expression'))),
    tabPanel("Plots",
             sidebarPanel(
                 downloadButton('downloadHeatmap', 'Download Heatmap'))),
    
    
    
    
    
    
    )
    
    
    
)
)
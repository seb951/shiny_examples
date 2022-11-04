#
library(shiny)
library(cluster)
library(DT)
library(wesanderson)
library(shinythemes)

shinyUI(navbarPage(title = "Clustering example",theme = shinytheme("cerulean"),
                   header = list(
    strong("sebastien.renaut@gmail.com, fall 2022"),
    htmlOutput("general")),
    tabPanel("Kmeans",
             sidebarPanel(
               fileInput("GEXP_file","Gene Expression dataset (.csv with header)", accept=".csv", width = NULL, placeholder = ""),
               fileInput("CLINICAL_file","Clinical dataset (.csv with header)", accept=".csv", width = NULL, placeholder = ""),
               sliderInput("k_selected","Select a K based on Silhouette plots",min = 2, max = 8, value = 2),
               checkboxGroupInput("pc", "Choose two PCs for visualisation",
                                  choiceNames = paste0("PC",1:8),
                                  selected = c(1,2),
                                  choiceValues =  1:8)
             ),
             mainPanel(h4("Kmeans"),
                       htmlOutput("kmeans"),
                       plotOutput("sil"),
                       plotOutput("detailedsil"),
                       plotOutput("pca"))
    ),
    
    tabPanel("Hierarchical",
             sidebarPanel(
               sliderInput("k_selected_hc","Select a K to cut tree",min = 2, max = 8, value = 2),
             ),
             mainPanel(
                 h4("Hierchical Clustering"),
               htmlOutput("dendrogram"),
               plotOutput("heatmap"),
               plotOutput("detailedsil_hc")
             )
    ),
    tabPanel("Summary",
             sidebarPanel(
               radioButtons("variable", "Choose a clinical variable to display:",
                            selected = "age",
                            choices= c("kmeans","dendrogram","age","stage","histology","gender","survival"))
             ),
                mainPanel(
                  h4("The Clinical Dataset"),
                 htmlOutput("clinical"),
                   plotOutput("summary_data")),
             mainPanel(
                 downloadLink('downloadData', 'Download Clinical Data'),
                 DT::dataTableOutput("mesotable"))
    )
    
    
    
)
)
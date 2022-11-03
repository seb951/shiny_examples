#
library(shiny)
library(cluster)
library(DT)
library(wesanderson)

shinyUI(fluidPage(
    
    # Application title
    titlePanel("Clustering example: mesothelioma dataset"),
    h4("sebastien.renaut@gmail.com"),
    htmlOutput("toto"),
    h4("fall 2022"),
    htmlOutput("general"),
    sidebarLayout(
        sidebarPanel(
            fileInput("GEXP_file","Gene Expression dataset (.csv with header)", accept=".csv", width = NULL, placeholder = ""),
            fileInput("CLINICAL_file","Clinical dataset (.csv with header)", accept=".csv", width = NULL, placeholder = "")
        ),
        mainPanel(
            textInput("dataset_name", "What is your dataset called? ", "mesothelioma")
        )
    ),
    
    
    
    h4("The Gene Expression Dataset"),
    htmlOutput("gene_expression"),
    
    
      sidebarLayout(
        sidebarPanel(
        #    sliderInput("krange","Select a K range to calculate Silhouette",min = 2, max = 8, value = c(2, 5))),
          sliderInput("k_selected","Select a K based on both Silhouette plots",min = 2, max = 8, value = 2)),
         mainPanel(
                plotOutput("sil")
            )),
    sidebarLayout(
      sidebarPanel(),
        mainPanel(
          plotOutput("detailedsil")
      )),
    sidebarLayout(
        sidebarPanel(
            checkboxGroupInput("pc", "Choose two PCs for visualisation",
                               choiceNames = paste0("PC",1:8),
                               selected = c(1,2),
                               choiceValues =  1:8
            ),
        ),
        mainPanel(
            plotOutput("pca")
        )),
    

    
    h4("The Clinical Dataset"),
    htmlOutput("clinical"),
    DT::dataTableOutput("mesotable"),
    
    h4("Summmarized data below:"),
    sidebarLayout(
        sidebarPanel(
            radioButtons("variable", "Choose a clinical variable to display:",selected = "age",choices= c("cluster","age","stage","histology","gender","survival"))),
        mainPanel(
            plotOutput("summary_data")
        )),
    
    
    
)
)
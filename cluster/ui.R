#
library(shiny)
library(cluster)
library(DT)
library(wesanderson)

shinyUI(fluidPage(
    
    # Application title
    titlePanel("Clustering example: mesothelioma dataset"),
    h4("sebastien.renaut@gmail.com"),
    h4("fall 2022"),
    htmlOutput("general"),
    
    h4("The Mesothelioma Gene Expression Dataset"),
    htmlOutput("gene_expression"),
    sidebarLayout(
        sidebarPanel(
            sliderInput("krange","Select a K range to calculate Silhouette",min = 2, max = 8, value = c(2, 5))),
            mainPanel(
                plotOutput("sil")
            )),
    
    sidebarLayout(
        sidebarPanel(
            sliderInput("k_selected","Select a K based on Silhouette plot above",min = 2, max = 8, value = 2),
            #radioButtons("pcX", "Choose two PCs:",
             #        choiceNames = paste0("PC",1:10),
              #       selected = c(1),
               #      choiceValues = 1:10
                #     )
            checkboxGroupInput("pc", "Choose two PCs",
                               choiceNames = paste0("PC",1:8),
                               selected = c(1,2),
                               choiceValues =  1:8
            ),
        ),
        mainPanel(
            plotOutput("pca")
        )),
    

    
    h4("The Mesothelioma Clinical Dataset"),
    htmlOutput("clinical"),
    DT::dataTableOutput("mesotable"),
    
    h4("Summmarized data below:"),
    sidebarLayout(
        sidebarPanel(
            radioButtons("variable", "Choose a clinical variable to display:",selected = "age",choices= c("cluster","age","stage","histology","gender","survival"))),
        mainPanel(
            plotOutput("summary_data")
        ))
    
    
)
)
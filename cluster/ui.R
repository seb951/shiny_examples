#
library(shiny)
library(cluster)
library(DT)


shinyUI(fluidPage(
    
    # Application title
    titlePanel("Clustering example: mesothelioma gene expression dataset"),
    h4("sebastien.renaut@gmail.com"),
    h4("fall 2022"),
    htmlOutput("text_cluster"),
    
    h4("The Mesothelioma Clinical Dataset"),
    DT::dataTableOutput("mesotable"),
    
    h4("Summmarized data below:"),
    sidebarLayout(
        sidebarPanel(
            radioButtons("variable", "Choose a clinical variable to display:",selected = "years_to_birth",choices= c("years_to_birth","pathologic_stage","histological_type","gender","overall_survival"))),
        mainPanel(
            plotOutput("summary_data")
        )),
    
    
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
    

    
)
)